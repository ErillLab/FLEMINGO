import random
import numpy as np
import math
import copy
global null_models
null_models = {}

def norm_cdf(x, mu, sigma):
    ''' Cumulative distribution function for the normal distribution. '''
    if sigma == 0:
        if x < mu:
            return 0
        elif x > mu:
            return 1
        elif x == mu:
            return 0.5
    else:
        z = (x-mu)/abs(sigma)  # z-score
        return (1.0 + math.erf(z / math.sqrt(2.0))) / 2.0

def norm_pf(x, y, mu, sigma):
    """
    Probablity function for normal distribution.
    Considering that the observed x is an integer, the probability of x is
    defined as the probability of observing a value within x - 0.5 and x + 0.5,
    given a normal distribution specified by the given mu and sigma.
    """
    if sigma != 0:
        p = norm_cdf(y, mu, sigma) - norm_cdf(x, mu, sigma)
    else:
        # when sigma is 0
        if x == mu:
            p = 1
        else:
            p = 0
    return p


class ShapeObject:

    def __init__(self, *args) -> None:
        if len(args) == 4:
            self.constructor_1(args[0], args[1], args[2], args[3])
        elif len(args) == 5:
            self.constructor_2(args[0], args[1], args[2], args[3], args[4])
        else:
            raise Exception("Wrong ShapeObject constructor arguments.")

    def constructor_1(self, rec_type, rec_size, mu = None, sigma = None):
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
            mu: mean for distribution
            sigma: standard deviation for distribution

        Returns:
            Fully operational shape object ready to be placed
        """
        self.type = rec_type
        self.length = rec_size
        
        # The following attributes are set by set_null_model()
        self.null_model = []
        self.edges = []
        self.min_mu = None
        self.max_mu = None
        self.set_null_model()  # Sets attributes above
        
        # Set mu and sigma
        self._mu = mu
        self._sigma = sigma
        
        if mu == None:
            self._set_random_mu()
        if sigma == None:
            self._set_random_sigma()
        
        ## config file params
        self.pseudo_count = 1e-100

        self.alt_model = []
        # array for alt model must be filled prior to placement
        self.set_alt_model()

    
    def constructor_2(self, rec_type, rec_size, config, mu = None, sigma = None):
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
        
        # The following attributes are set by set_null_model()
        self.null_model = []
        self.edges = []
        self.min_mu = None
        self.max_mu = None
        self.set_null_model()  # Sets attributes above
        
        # Set mu and sigma
        # self.min_mu = self.edges[0]  # !!! redundant code (happens in set_null_model)
        # self.max_mu = self.edges[-1]  # !!! redundant code (happens in set_null_model)
        self._mu = mu
        self._sigma = sigma
        
        if mu == None:
            self._set_random_mu()
        if sigma == None:
            self._set_random_sigma()
        
        self.pseudo_count = config["PSEUDO_COUNT"]
        self.probability_mu_mutation = config["PROBABILITY_MU_MUTATION"]
        self.probability_sigma_mutation = config["PROBABILITY_SIGMA_MUTATION"]
        self.magnitude_mu_mutation = config["MAGNITUDE_MU_MUTATION"]
        self.magnitude_sigma_mutation = config["MAGNITUDE_SIGMA_MUTATION"]
        self.mu_mutator = config["MU_MUTATOR"]
        self.sigma_mutator = config["SIGMA_MUTATOR"]
        self.probability_change_rec_type = config["PROBABILITY_CHANGE_REC_TYPE"]
        self.probability_increase_length = config["PROBABILITY_INCREASE_LENGTH"]
        self.probability_decrease_length = config["PROBABILITY_DECREASE_LENGTH"]
        self.min_length = config["MIN_LENGTH"]
        self.max_length = config["MAX_LENGTH"]
        self.num_bins = config["NUM_BINS"]
        
        self.alt_model = []
        # array for alt model must be filled prior to placement
        self.set_alt_model()
    
    def _set_random_mu(self):
        '''
        Chooses a random value for mu during Shape object initialization. The
        value is chosen uniformly from the possible values for the given shape.
        This function is called by the __init__ function.
        '''
        self._mu = np.random.uniform(self.min_mu, self.max_mu)
    
    def _set_random_sigma(self):
        '''
        Chooses a random value for sigma during Shape object initialization.
        Sigma can be any non-negative real number. The value is drawn from an
        exponential distribution. The parameter of the exponential is chosen so
        that the expected value of sigma (avg_sigma) is such that the interval
        from mu-3*avg_sigma to mu+3*avg_sigma spans as much as half of the total
        range. This means that the interval from mu to mu+3*avg_sigma must span
        one quarter of the total range:
            range/4 = mu+3*avg_sigma - mu
            range/4 = 3*avg_sigma
        Therefore,
            avg_sigma = range/12
        
        This function is called by the __init__ function.
        '''        
        avg_sigma = (self.edges[-1] - self.edges[0]) / 12
        # Random Generator
        rng = np.random.default_rng()
        # Draw value of sigma from exponential distribution
        self._sigma = rng.exponential(avg_sigma)
    
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
        self.edges = null_models[self.type][self.length]["bins"]
        self.min_mu = self.edges[0]
        self.max_mu = self.edges[-1]

    def set_alt_model(self):
        """Sets the array of frequencies corresponding to the distribution of the
        shape recognizer 
        """
        
        # Make sure mu is within bounds (min and max observed values in the null)
        if self._mu < self.min_mu:
            self._mu = self.min_mu
        elif self._mu > self.max_mu:
            self._mu = self.max_mu
        
        self.alt_model = []

        # bin array is one longer than frequency array
        # computed in the same way that null model for connectors is computed

        auc = norm_cdf(self.edges[-1], self._mu, self._sigma) - norm_cdf(self.edges[0], self._mu, self._sigma)

        auc += (len(self.edges) - 1) * self.pseudo_count

        log_auc = np.log(auc)
        
        for i in range(0, len(self.edges) - 1):
            score = norm_pf(self.edges[i], self.edges[i + 1], self._mu, self._sigma)
            self.alt_model.append(np.log(score + self.pseudo_count) - log_auc)

    def mutate(self, org_factory):
        '''
        Mutation for a Shape recognizer. The two parameters (mu and sigma) can be
        mutated (depending on the respective mutation rates).

        Args:
            org_factory: Organism Facory
            ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            !!! Unused argument (just for consistency with the pssm mutate function)
            ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        '''
        # keeps track of if a mutation happend. If it did, alternative model is updated
        mutated = False
        
        # Code to keep track if mutations that shift the boundaries of the
        # recognizer placement have occurred
        displacement_code = [0, 0]  # left and right displacement (in bp)
        
        # Mutate sigma
        if random.random() < self.probability_sigma_mutation:
            self._mutate_sigma()
            mutated = True
        
        # Mutate mu
        if random.random() < self.probability_mu_mutation:
            self._mutate_mu()
            mutated = True

        # size mutations, bounded by min and max lengths specified by config
        # null models must be updated whenever size changes
        if random.random() < self.probability_increase_length and self.length < self.max_length:
            self._increase_length(displacement_code)
            mutated = True

        if random.random() < self.probability_decrease_length and self.length > self.min_length:
            self._decrease_length(displacement_code)
            mutated = True
        
        if mutated:
            self.set_alt_model()
        
        # If the recognizer boundaries have changed, report the `displacement_code`
        # so that connectors can eventually be adjusted if necessary
        if displacement_code != [0, 0]:
            return displacement_code
    
    def _mutate_sigma(self):
        '''
        Change the value of the sigma parameter. The way it's modified depends
        on the `sigma_mutator` parameter specified in the config file.
        
        =======
        Warning
        =======
        This function is meant to be called by the `mutate` function.
        If you call this function outside the `mutate` function, the values
        for each bin of the 'alternative model' will be incorrect unless you
        update them by calling:
            
            self.set_alt_model()
        
        '''
        if self.sigma_mutator=="linear":
            # Update sigma with a random permutation within allowed interval
            self._sigma = abs(
                self._sigma + random.uniform(-self.magnitude_sigma_mutation,
                                             self.magnitude_sigma_mutation)
            )
        elif self.sigma_mutator=="log":
            
            # !!! Special case with hard-coded pseudosigma
            if self._sigma == 0:
                self._sigma = 0.00001
            else:
                base = self.magnitude_sigma_mutation
                logb_sigma = np.log(self._sigma) / np.log(base)
                shift = random.uniform(-1, 1)
                # Apply a shift in the range (-1, 1) to the log-sigma
                logb_sigma += shift
                self._sigma = base**logb_sigma
    
    def _mutate_mu(self):
        '''
        Change the value of the mu parameter. The way it's modified depends
        on the `mu_mutator` parameter specified in the config file.
        
        =======
        Warning
        =======
        This function is meant to be called by the `mutate` function.
        If you call this function outside the `mutate` function, the values
        for each bin of the 'alternative model' will be incorrect unless you
        update them by calling:
            
            self.set_alt_model()
        
        '''
        if self.mu_mutator=="linear":
            # Update mu with a random permutation within allowed interval
            self._mu = self._mu + random.uniform(-self.magnitude_mu_mutation, self.magnitude_mu_mutation)

        elif self.mu_mutator=="log":
            base = self.magnitude_mu_mutation
            logb_mu = np.log(self._mu) / np.log(base)
            shift = random.uniform(-1, 1)
            # Apply a shift in the range (-1, 1) to the log-mu
            logb_mu += shift
            self._mu = base**logb_mu
        
        elif self.mu_mutator=="standard":
            self._mu = random.gauss(self._mu, self._sigma)
    
    def _increase_length(self, displacement_code, token_arg=None):
        '''
        Increase length of shape recognizer.
        The `displacement_code` is updated to report the shift of the bounds of
        the recognizer. This info is used to adjust the adjacent connectors
        accordingly.
        Example: when the recognizer grows by one bp, one of the adjacent
        connectors may have to reduce its `mu` by one (it may not happen if the
        shape recognizer is a terminal node).
        
        =======
        Warning
        =======
        This function is meant to be called by the `mutate` function.
        If you call this function outside the `mutate` function, the values
        for each bin of the 'alternative model' will be incorrect unless you
        update them by calling:
            
            self.set_alt_model()
        
        '''
        # do only if allowed
        # Checking this condition is useful when calling this function from
        # outside the `mutate` function
        if self.length < self.max_length:
            self.length += 1 
            self.set_null_model()
            
            # Compile the `displacement_code`. This allows to adjust the mu of
            # adjacent connectors if necessary. Differently from the PSSM, which
            # grows by adding a column to the left or to the right, the growth
            # of the shape recognizer is not explicitly directional. Therefore,
            # we choose randomly which connector (right or left) will be asked
            # to shrink `mu` by 1 unit.
            if random.random() < 0.5:
                # The left bound ([0]) of the recognizer shifts 1 bp to the left (-1)
                displacement_code[0] -= 1
            else:
                # The right bound ([1]) of the recognizer shifts 1 bp to the right (+1)
                displacement_code[1] += 1
    
    def _decrease_length(self, displacement_code):
        '''
        Decrease length of shape recognizer.
        The `displacement_code` is updated to report the shift of the bounds of
        the recognizer. This info is used to adjust the adjacent connectors
        accordingly.
        Example: when the recognizer shrinks by one bp, one of the adjacent
        connectors may have to increase its `mu` by one (it may not happen if
        the shape recognizer is a terminal node).
        
        =======
        Warning
        =======
        This function is meant to be called by the `mutate` function.
        If you call this function outside the `mutate` function, the values
        for each bin of the 'alternative model' will be incorrect unless you
        update them by calling:
            
            self.set_alt_model()
        
        '''
        # do only if allowed
        # Checking this condition is useful when calling this function from
        # outside the `mutate` function
        if self.length > self.min_length:
            self.length -= 1 
            self.set_null_model()
            
            # Compile the `displacement_code`. This allows to adjust the mu of
            # adjacent connectors if necessary. Differently from the PSSM, which
            # shrinks by dropping a column to the left or to the right, the
            # reduction of the shape recognizer is not explicitly directional.
            # Therefore, we choose randomly which connector (right or left) will
            # be asked to increase `mu` by 1 unit.
            if random.random() < 0.5:
                # The left bound ([0]) of the recognizer shifts 1 bp to the right (+1)
                displacement_code[0] += 1
            else:
                # The right bound ([1]) of the recognizer shifts 1 bp to the lefet (-1)
                displacement_code[1] -= 1
    
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
    
    def is_connector(self) -> bool:
        return False
    
    def is_pssm(self) -> bool:
        return False
    
    def is_shape(self) -> bool:
        return True

    def get_type(self):
        """Returns the one character definiton of a recognizers type
        used for generating string of recognizer types.
        """
        if self.type == "mgw":
            return 'm'
        if self.type == "prot":
            return 't'
        if self.type == "helt":
            return 'h'
        if self.type == "roll":
            return 'r'

    def to_json(self) -> dict:
        """Get Shape object in JSON format

        Returns:
            shape in JSON format
        """
        shape = {}
        shape["objectType"] = "shape"
        shape["recType"] = self.type
        shape["mu"] = self._mu
        shape["sigma"] = self._sigma
        shape["length"] = self.length
        return shape
