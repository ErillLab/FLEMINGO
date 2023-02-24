import random
import numpy as np
import math
import copy

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
    def __init__(self, rec_type, rec_size, mu, sigma, config, null_models):
        self.type = rec_type
        self.length = rec_size
        self._mu = mu
        self._sigma = sigma
        self.null_model = []
        self.bins = []
        self.alt_model = []
        self.nulls = null_models
        self.set_null_model(null_models)
        self.set_alt_model()

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

    def set_null_model(self, null_models):
        rec_type = ''
        if self.type == 'm':
            rec_type = "mgw"
        if self.type == 't':
            rec_type = "prot"
        if self.type == 'r':
            rec_type = "roll"
        if self.type == 'h':
            rec_type = "helt"
        self.null_model = null_models[rec_type][self.length]["values"]
        self.bins = null_models[rec_type][self.length]["bins"]

    def set_alt_model(self):
        alt_model = []
        for i in self.bins:
            score = norm_pf(i + 0.05, self._mu, self._sigma)
            if score < 0.001:
                score = 0.0001
            alt_model.append(score)
        self.alt_model = np.array(alt_model, dtype=np.dtype("f"))

    def mutate(self, org_fac):
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
                
       
        # MU MUTATION
        if random.random() < self.mutate_probability_mu:
            #determine type of mutation (linear or log)
            if self.mu_mutator=="linear":
                # Update mu with a random permutation within allowed interval
                self._mu = abs(
                    self._mu + random.uniform(-self.mutate_variance_mu,
                                              self.mutate_variance_mu)
                )
            elif self.mu_mutator=="log":
                base = self.mutate_variance_mu
                logb_mu = np.log(self._mu) / np.log(base)
                shift = random.uniform(-1, 1)
                # Apply a shift in the range (-1, 1) to the log-mu
                logb_mu += shift
                self._mu = base**logb_mu
            
            elif self.mu_mutator=="standard":
                self._mu = abs(random.gauss(self._mu, self._sigma))

        if random.random() < self.mutate_probability_increase_size and self.length < self.max_columns:
            self.length += 1 
            self.set_null_model(self.nulls)
            self.set_alt_model()

        if random.random() < self.mutate_probability_decrease_size and self.length > self.min_columns:
            self.length -= 1
            self.set_null_model(self.nulls)
            self.set_alt_model()

    def export(self, export_file) -> None:
        """Exports pssm to a file

        Args:
            export_file: File to write the output
        """
        
        
        export_file.write("\n" + str(self._mu) + " " + str(self._sigma) + " " + self.type)


