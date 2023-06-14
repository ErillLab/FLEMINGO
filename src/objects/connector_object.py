"""C object
   Connector objects connect recognizer objects using two parameters:
   - mu: determines the ideal distance between the two recognizer objects
   - sigma: determines the standard deviation if the distance between the two
     recognizer objects
   The "energy" contribution of a connector object depends on the
   relative distance of the two recognizer nodes that it connects when
   placed on the sequence.
"""

import random
import numpy as np
import math


def g(x, mu, sigma):
    return -(1/2) * ((x-mu)/sigma)**2

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

def prob_of_d(d, L, N):
    ''' Given N randomly chosen integers in [1,L], this function returns the
    probability that two consecutive integers (after sorting) are found at a
    distance d. '''
    if 1 <= d and d <= L-N+1:
    #if 1 <= d and d <= L:
        num = math.comb(L-d, N-1)
        den = math.comb(L, N)
        return num / den
    else:
        # !!! This probability should be 0. In the python implementation of the
        # placement algorithm "impossible" paths are considered: gaps that are
        # too large (some PSSM would not fit the DNA sequence on the right) are
        # not skipped! Instead the values are stored in the matrix. As in the
        # Needleman–Wunsch algorithm, they just will never reach the bottom line
        # of the alignment/placement matrix, so they are implicitly discarded.
        # No placement with such large gap can ever be observed in any practical
        # situations (because PSSMs are not allowed to overlap, no matter what).
        # However, we can't return 0, because it would cause a zero division
        # error when computing those useless cells in the upper right corner of
        # the placement matrix. Here we return 10^-10 instead, so that no error
        # is encountered when filling the matrix.
        return 10**(-10)


class ConnectorObject():
    """Connector Object is a node that connects two recognizer objects
    """
    
    def __init__(self, _mu: int, _sigma: int, config: dict, max_seq_length: int):
        """Connector constructor: 
            - gets/sets mu and sigma
            - sets all connector-specific configuration items

        Args:
            _mu: Mean distance between node1 and node2
            _sigma: Variance of distance between node 1 and node2
            config: Node-level configuration specs as loaded from config.json

        """
        # set mu and sigma
        self._mu = _mu  # Mean discance between the connected nodes
        self._sigma = _sigma  # Standard deviation of the distance
        # set connector-specific configuration parameters
        self.mutate_probability_sigma = config["MUTATE_PROBABILITY_SIGMA"]
        self.mutate_probability_mu = config["MUTATE_PROBABILITY_MU"]
        self.mutate_variance_sigma = config["MUTATE_VARIANCE_SIGMA"]
        self.mutate_variance_mu = config["MUTATE_VARIANCE_MU"]
        self.max_seq_length = max_seq_length
        self.sigma_mutator = config["SIGMA_MUTATOR"] #log or linear
        self.mu_mutator = config["MU_MUTATOR"] #log or linear
        self.pseudo_count = config["PSEUDO_COUNT"]
        
        # precompute connector energies for expected length range
        self.stored_pdfs = []
        self.stored_cdfs = []
        self.set_precomputed_pdfs_cdfs()

        self.adjust_score_threshold = 0
    
    # Setters
    def set_mu(self, _mu: int) -> None:
        """Set mu variable

        Args:
            _mu: Mean distance between nodes connected by connector
        """
        self._mu = _mu

    def set_sigma(self, sigma: int) -> None:
        """Set sigma variable

        Args:
            sigma: Standard deviation of the distance between nodes connected
            by connector
        """
        self._sigma = sigma
    
    def set_precomputed_pdfs_cdfs(self) -> None:
        """Set stored_pdfs variable and stored_cdfs variable. This is done
        through computation of cdfs for all intervals with width 1 from -0.5 to 
        max_seq_length + 0.5. the pdfs can also be computed using these
        computations by taking the difference a cdf and its previous cdf.
        a pseudo count is added to each interval in the pdfs. The cdfs have
        a pseudo count of (pseudo_count * j) added where j is the number of
        pdf bins that have had the pseudo count added so far.
        """
        
        # Delete previous values
        self.stored_pdfs = [None] * self.max_seq_length
        self.stored_cdfs = [None] * self.max_seq_length
       
        prev_cdf = norm_cdf(-0.5, self._mu, self._sigma)
        cdf_0 = prev_cdf
        offset = 0
        # range is 1 to max_seq_length + 1 because first pdf is cdf(0.5) - cdf(-0.5)
        # and the last pf should be cdf(max_seq_length + 0.5) - cdf(max_seq_length - 0.5)
        # since we're starting at 1, for each index in the auc array we need
        # to actually use the previous iteration of the auc since the first index
        # should be cdf(-0.5) - cdf(-0.5)

        for j in range(1, self.max_seq_length + 1):
            cdf = norm_cdf(j - 0.5, self._mu, self._sigma)
            auc = np.log2(prev_cdf - cdf_0 + (self.pseudo_count * j))
            pf = np.log2(cdf - prev_cdf + self.pseudo_count)

            self.stored_pdfs[offset] = np.double(pf)
            self.stored_cdfs[offset] = np.double(auc)
            offset += 1
            prev_cdf = cdf
    
    def mutate(self, org_factory) -> None:
        """mutation for a connector

        Args:
            org_factory(organism_factory): Organism Facory
        """
        
        mutated = False
        
        # Mutate sigma
        if random.random() < self.mutate_probability_sigma:
            self._mutate_sigma()
            mutated = True
            
            # #determine type of mutation (linear or log)
            # if self.sigma_mutator=="linear":
            #     # Update sigma with a random permutation within allowed interval
            #     self._sigma = abs(
            #         self._sigma + random.uniform(-self.mutate_variance_sigma,
            #                                      self.mutate_variance_sigma)
            #     )
            # elif self.sigma_mutator=="log":
            #     base = self.mutate_variance_sigma
            #     logb_sigma = np.log(self._sigma) / np.log(base)
            #     shift = random.uniform(-1, 1)
            #     # Apply a shift in the range (-1, 1) to the log-sigma
            #     logb_sigma += shift
            #     self._sigma = base**logb_sigma
        
        # Mutate mu
        if random.random() < self.mutate_probability_mu:
            self._mutate_mu()
            mutated = True
            
            # #determine type of mutation (linear or log)
            # if self.mu_mutator=="linear":
            #     # Update mu with a random permutation within allowed interval
            #     self._mu = abs(
            #         self._mu + random.uniform(-self.mutate_variance_mu,
            #                                   self.mutate_variance_mu)
            #     )
            # elif self.mu_mutator=="log":
            #     base = self.mutate_variance_mu
            #     logb_mu = np.log(self._mu) / np.log(base)
            #     shift = random.uniform(-1, 1)
            #     # Apply a shift in the range (-1, 1) to the log-mu
            #     logb_mu += shift
            #     self._mu = base**logb_mu
            
            # elif self.mu_mutator=="standard":
            #     self._mu = abs(random.gauss(self._mu, self._sigma))
        
        if mutated:
            # Recompute PDF and CDF values
            self.set_precomputed_pdfs_cdfs()
    
    def _mutate_sigma(self):
        '''
        Change the value of the sigma parameter. The way it's modified depends
        on the `sigma_mutator` parameter specified in the config file.
        
        =======
        Warning
        =======
        This function is meant to be called by the `mutate` function.
        If you call this function outside the `mutate` function, the pdfs and
        cdfs values will be incorrect unless you update them by calling:
            
            self.set_precomputed_pdfs_cdfs()
        
        '''
        if self.sigma_mutator=="linear":
            # Update sigma with a random mutation within allowed interval
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
    
    def _mutate_mu(self):
        '''
        Change the value of the mu parameter. The way it's modified depends
        on the `mu_mutator` parameter specified in the config file.
        
        =======
        Warning
        =======
        This function is meant to be called by the `mutate` function.
        If you call this function outside the `mutate` function, the pdfs and
        cdfs values will be incorrect unless you update them by calling:
            
            self.set_precomputed_pdfs_cdfs()
        
        '''
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
    
    def null_gap_likelihood(self, gap_size, recog_sizes, seq_len):
        # For each recog of size s, we must subtract (s-1). E.g., for a recog of
        # size 5 we must subtract 4.
        effective_len = seq_len - sum(recog_sizes) + len(recog_sizes)
        
        # Distance = gap + 1
        # The length to be used as input is the effective length
        # Number of recognizers is the length of the list of recog sizes
        return prob_of_d(gap_size+1, effective_len, len(recog_sizes))

    def get_numerator(self, d, s_dna_len, recog_sizes) -> float:
        if d<self.max_seq_length:
            numerator = self.stored_pdfs[d]
        else:
            numerator = norm_pf(d, self._mu, self._sigma)
        
        # Normalize by AUC within the range of observable d values
        max_d = s_dna_len - 1  # Maximum d observable
        # min_d = -1 * max_d  # Minimum d observable
        
        if self._sigma == 0:
            auc = 1.0  # all the gaussian is within the (-(L-1), +(L-1)) range
        else:
            if max_d<self.max_seq_length:
                auc = self.stored_cdfs[max_d] - self.stored_cdfs[0]
            else:
                auc = (norm_cdf(max_d, self._mu, self._sigma) -
                       norm_cdf(0, self._mu, self._sigma))
        
        # Avoid zero-division error
        # This will never happen, unless an organism evolves an extremely large sigma
        if auc < 1e-100:
            auc = 1e-100 
            print("AUC was 0 with mu =", self._mu, "and sigma =", self._sigma)                
        
        # avoid log(0) error when computing e_connector
        if numerator < 1e-100:
            numerator = 1e-100        

        # Apply normalization
        return np.log2(numerator / auc)

    def get_score(self, d, s_dna_len, recog_sizes) -> float:
        """ Returns the score of the connector, given the observed distance
            and the length of the DNA sequence on which it is being evaluated.
            Parameters
            ----------
            d : observed distancce
            s_dna_len : length of DNA sequence on which connector is placed
    
            Returns
            -------
            e_connector : the connector's score
            
            The score of the connector is computed as a log-likelihood ratio.
            The numerator is the probability of observing the distance provided
            given the connector's parameters.
            The probability of observing distance d, given the connector's
            parameters is given by norm_pf.
            The probability is then normalized by the cumulative probability
            function within the observable range on the sequence.
            
            The denominator is the probability of observing the distance
            provided under the null hypothesis.
            
            To speed up the process, pdfs and cdfs are precomputed for the
            "expected" distance range during connector construction, and looked
            up, rather than recomputed here.

        """
        # Numerator
        if d<self.max_seq_length:
            numerator = self.stored_pdfs[d]
        else:
            numerator = norm_pf(d, self._mu, self._sigma)
        
        # Normalize by AUC within the range of observable d values
        max_d = s_dna_len - 1  # Maximum d observable
        # min_d = -1 * max_d  # Minimum d observable
        
        if self._sigma == 0:
            auc = 1.0  # all the gaussian is within the (-(L-1), +(L-1)) range
        else:
            if max_d<self.max_seq_length:
                auc = self.stored_cdfs[max_d] - self.stored_cdfs[0]
            else:
                auc = (norm_cdf(max_d, self._mu, self._sigma) -
                       norm_cdf(0, self._mu, self._sigma))
        
        # Avoid zero-division error
        # This will never happen, unless an organism evolves an extremely large sigma
        if auc < 1e-100:
            auc = 1e-100 
            print("AUC was 0 with mu =", self._mu, "and sigma =", self._sigma)                
        
        # avoid log(0) error when computing e_connector
        if numerator < 1e-100:
            numerator = 1e-100        

        # Apply normalization
        numerator = numerator / auc
        
        # Denominator
        # The denominator is p(d) according to the null model
        denominator = self.null_gap_likelihood(abs(d), recog_sizes, s_dna_len)
        
        # compute additive connector energy term as log-likelihood ratio
        e_connector = np.log2(numerator / denominator)
        
        return e_connector

    def print(self, debug=False) -> None:
        """Prints the connector mu and sigma values
        """
        print(" m: {} s: {}".format(self._mu, self._sigma))
        if debug:
            print(" adjust_score_threshold: {}".format(self.adjust_score_threshold))

    def export(self, export_file) -> None:
        """Exports Connector data to the given file

        Args:
            export_file (file): File to export the conector

        """
        export_file.write("\n m: {} s: {}".format(self._mu, self._sigma))

    def is_connector(self) -> bool:
        """node is connector

        Returns:
            True because is a connector
        """
        return True

    def is_pssm(self) -> bool:
        """node is not a pssm

        Returns:
            False because is a connector
        """
        return False

