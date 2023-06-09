import random
import numpy as np
import math
import copy
import models.models as models
global null_models
null_models = {}

def norm_cdf(x, mu, sigma):
    ''' Cumulative distribution function for the normal distribution. '''
    z = (x-mu)/abs(sigma)
    return (1.0 + math.erf(z / math.sqrt(2.0))) / 2.0

def norm_pf(x, mu, sigma):
    """
    Probablity function for normal distribution.
    Considering that the observed x is an integer, the probability of x is
    defined as the probability of observing a value within x - 0.5 and x + 0.5,
    given a normal distribution specified by the given mu and sigma.
    """
    if sigma != 0:
        p = norm_cdf(x+0.5, mu, sigma) - norm_cdf(x-0.5, mu, sigma)
    else:
        # when sigma is 0
        if x == mu:
            p = 1
        else:
            p = 0
    return p


class ShapeObject:
    def __init__(self, rec_type, rec_size, config, mu = None, sigma = None):
        """ShapeObject: recognizers a specific DNA-shape {mgw, prot, roll, 
        or helt}

        All null models ever computed are stored in the global null_models 
        variable. Each shape has exactly one null model corresponding to 
        the DNA-shape feature that it recognizers and its length.

        Each shape has exactly one alternative model represented by a 
        normal distribution following a mean, mu, and standard deviation, sigma. 
        It's representation distribution is calculated of the same intervals 
        that the null models are.

        Shape features can be mutated by either increasing or decreasing its mu, 
        its sigma, and/or its length. mu is bounded by the min and max obseved
        scores in the null model, sigma is any non-negative number, and length is
        bounded by the config file.

        Note: minimum length can be modified but can NOT be less than 5.

        Args:
            rec_type: kind of shape feature
            rec_size: length of the recognizer
            config: shape recognizer config parameters
            mu: mean for distribution
            sigma: standard deviation for distribution

        Returns:
            Fully operational shape object ready to be placed
        """
        self.type = rec_type
        self.length = rec_size

        self.null_model = []
        self.bins = []
        self.alt_model = []
        self.set_null_model()

        self.min_mu = self.bins[0]
        self.max_mu = self.bins[-1]
        self._mu = mu
        self._sigma = sigma

        self.pseudo_count = config["PSEUDO_COUNT"]
        self.mutate_probability_sigma = config["MUTATE_PROBABILITY_SIGMA"]
        self.mutate_probability_mu = config["MUTATE_PROBABILITY_MU"]
        self.sigma_mutator = config["SIGMA_MUTATOR"]
        self.mu_mutator = config["MU_MUTATOR"]
        self.mutate_variance_sigma = config["MUTATE_VARIANCE_SIGMA"]
        self.mutate_variance_mu = config["MUTATE_VARIANCE_MU"] 
        self.mutate_probability_change_rec_type = config["MUTATE_PROBABILITY_CHANGE_REC_TYPE"]
        self.mutate_probability_increase_size = config["MUTATE_PROBABILITY_INCREASE_SIZE"]
        self.mutate_probability_decrease_size = config["MUTATE_PROBABILITY_DECREASE_SIZE"]
        self.min_columns = config["MIN_COLUMNS"]
        self.max_columns = config["MAX_COLUMNS"]

        if mu == None or sigma == None:
            self._mu = np.random.uniform(self.min_mu, self.max_mu)
            self._sigma = np.random.uniform(0, 1)

        # array for alt model must be filled prior to placement
        self.set_alt_model()

        

    def set_null_model(self):
        """Sets the null model corresponding to shape recognizer's feature
        and length by retrieving information from the appropriate 
        sub-dictionary of null_models.

        Args:
            None

        Returns:
            None
        """
        self.null_model = null_models[self.type][self.length]["frequencies"]
        self.bins = null_models[self.type][self.length]["bins"]
        self.min_mu = self.bins[0]
        self.max_mu = self.bins[-1]

    def set_alt_model(self):
        """Sets the array of frequencies corresponding to the distribution of the
        shape recognizer 
        """

        self.alt_model = []

        # bin array is one longer than frequency array
        # computed in the same way that null model for connectors is computed
        for i in range(0, len(self.bins) - 1):
            score = norm_pf(i + 0.05, self._mu, self._sigma)
            self.alt_model.append(np.log(score + self.pseudo_count))

    def mutate(self, org_fac):
        """randomly mutates various attributes of the shape recognizer object.
        The recognizer may or may not be mutated

        Args:
            org_fac (has config information for determining mutations)

        Returns:
            None
            
        """
        # keeps track of if a mutation happend, if it did, then alternative
        # model is updated
        mutated = False

        # SIGMA MUTATION
        if random.random() < self.mutate_probability_sigma:
            #determine type of mutation (linear or log)
            if self.sigma_mutator=="linear":
                # Update sigma with a random permutation within allowed interval
                self._sigma = abs(
                    self._sigma + random.uniform(-self.mutate_variance_sigma,
                                                 self.mutate_variance_sigma)
                )
            elif self.sigma_mutator=="log":
                base = self.mutate_variance_sigma
                logb_sigma = np.log(self._sigma) / np.log(base)
                shift = random.uniform(-1, 1)
                # Apply a shift in the range (-1, 1) to the log-sigma
                logb_sigma += shift
                self._sigma = base**logb_sigma

            mutated = True
                
       
        # MU MUTATION
        if random.random() < self.mutate_probability_mu:
            #determine type of mutation (linear or log)
            if self.mu_mutator=="linear":
                # Update mu with a random permutation within allowed interval
                self._mu = self._mu + random.uniform(-self.mutate_variance_mu, self.mutate_variance_mu)

            elif self.mu_mutator=="log":
                base = self.mutate_variance_mu
                logb_mu = np.log(self._mu) / np.log(base)
                shift = random.uniform(-1, 1)
                # Apply a shift in the range (-1, 1) to the log-mu
                logb_mu += shift
                self._mu = base**logb_mu
            
            elif self.mu_mutator=="standard":
                self._mu = random.gauss(self._mu, self._sigma)

            mutated = True

        # size mutations, bounded by min and max lengths specified by config
        # null models must be updated whenever size changes
        if random.random() < self.mutate_probability_increase_size and self.length < self.max_columns:
            self.length += 1 
            self.set_null_model()
            mutated = True

        if random.random() < self.mutate_probability_decrease_size and self.length > self.min_columns:
            self.length -= 1
            self.set_null_model()
            mutated = True

        # mu is bounded by the minimum and maximum observed values in the
        # null
        if self._mu < self.min_mu:
            self._mu = self.min_mu

        if self._mu > self.max_mu:
            self._mu = self.max_mu

        if mutated:
            self.set_alt_model()


    def print(self) -> None:
        """Displays information about the shape recongizer

        Args:
            None

        Returns:
            None
            
        """
        print("Shape recognizer: ")
        print("Type:", self.type, "Length:", self.length)
        print("mu:", self._mu, "sigma", self._sigma)

    def export(self, export_file) -> None:
        """Exports shape to a file

        Args:
            export_file: File to write the output
        """
        export_file.write("\n" + str(self._mu) + " " + str(self._sigma) + " " + self.type + " " + str(self.length))

    def get_type(self):
        """Returns the one character definiton of a recognizers type
        used for generating string of recognizer types

        Args: 
            None

        Returns:
            None
        """
        if self.type == "mgw":
            return 'm'
        if self.type == "prot":
            return 't'
        if self.type == "helt":
            return 'h'
        if self.type == "roll":
            return 'r'
