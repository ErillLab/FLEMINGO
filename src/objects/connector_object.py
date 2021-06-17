"""C object
   Connector objects connect recognizer objects using two parameters:
   - mu: determines the ideal distance between the two recognizer objects
   - sigma: determines the standard deviation if the distance between the two
     recognizer objects
   The "energy" contribution of a connector object depends on the
   relative distance of the two recognizer nodes that it connects when
   placed on the sequence.
"""

# pylint: disable=E0402
# type: ignore
import random
import numpy as np
import math


def norm_pdf(x, mu, sigma):
    """ Normal probablity densitity function
        computes the normal probability density function
        used to compute the energy provided by the connector
        the distance at which connecting nodes are placed
        is the 'x' for the function, which will return the
        probability of observing said distance under the 
        connector model
   """
    if sigma != 0:
        var = float(sigma)**2
        denom = (2*math.pi*var)**.5
        num = math.exp(-(float(x)-float(mu))**2/(2*var))
        p = num/denom
    else:
        # when sigma is 0
        if x == mu:
            p = 1
        else:
            p = 0
    return p

def norm_cdf(x, mu, sigma):
    # Cumulative distribution function for the normal distribution
    z = (x-mu)/abs(sigma)
    return (1.0 + math.erf(z / math.sqrt(2.0))) / 2.0



# pylint: disable=R0902
class ConnectorObject():
    """Connector Object is a node that connects two recognizer objects
    """

    # pylint: disable=R0913
    def __init__(
            self,
            _mu: int,
            _sigma: int,
            config: dict,
    ):
        """Connector constructor: 
            - gets/sets mu and sigma
            - sets all connector-specific configuration items

        Args:
            _mu: Mean distance between node1 and node2
            _sigma: Variance of distance between node 1 and node2
            config: Node-level configuration specs as loaded from config.json

        """
        # set mu and sigma
        self._mu = _mu  # Mean discance
        self._sigma = _sigma  # Variance between elements
        # set connector-specific configuration parameters
        self.mutate_probability_sigma = config["MUTATE_PROBABILITY_SIGMA"]
        self.mutate_probability_mu = config["MUTATE_PROBABILITY_MU"]
        self.mutate_variance_sigma = config["MUTATE_VARIANCE_SIGMA"]
        self.mutate_variance_mu = config["MUTATE_VARIANCE_MU"]
        self.expected_seq_length = config["EXPECTED_SEQ_LENGTH"]
        self.sigma_mutator = config["SIGMA_MUTATOR"] #log or linear
        self.mu_mutator = config["MU_MUTATOR"] #log or linear
        
        # precompute connector energies for expected length range
        self.stored_pdfs = []
        self.stored_cdfs = []
        self.set_precomputed_pdfs_cdfs()
    
    # pylint: enable=R0913
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
        """Set stored_pdfs variable and stored_cdfs variable
        """
        
        # Delete previous values
        self.stored_pdfs = []
        self.stored_cdfs = []
        
        # Compute new values
        for dist in range(self.expected_seq_length):
            
            # Precompute PDF
            self.stored_pdfs.append(norm_pdf(dist, self._mu, self._sigma))
            
            # Precompute CDF
            if self._sigma != 0:
                self.stored_cdfs.append(norm_cdf(dist, self._mu, self._sigma))
            else:
                if dist<self._mu:
                    self.stored_cdfs.append(0.0)
                else:
                    self.stored_cdfs.append(1.0)
    
    # pylint: disable=W0613
    def mutate(self, org_factory) -> None:
        """mutation for a connector

        Args:
            org_factory(organism_factory): Organism Facory
        """
        
        # SIGMA MUTATION
        if random.random() < self.mutate_probability_sigma:
            #determine type of mutation (linear or log)
            if self.sigma_mutator=="linear":
                # Update sigma with a random permutation within allowed interval
                self._sigma = abs(
                    self._sigma + random.uniform(-self.mutate_variance_sigma,
                                                 self.mutate_variance_sigma)
                )
            else:
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
            else:
                base = self.mutate_variance_mu
                logb_mu = np.log(self._mu) / np.log(base)
                shift = random.uniform(-1, 1)
                # Apply a shift in the range (-1, 1) to the log-mu
                logb_mu += shift
                self._mu = base**logb_mu
        
        # Recompute PDF and CDF values
        self.set_precomputed_pdfs_cdfs()

    # pylint: enable=W0613

       
    def get_score(self, d, s_dna_len) -> float:
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
            parameters is given by norm_pdf.
            The probability is then normalized by the cumulative probability
            function within the observable range on the sequence.
            
            The denominator is the probability of observing the distance
            provided under the null hypothesis.
            
            To speed up the process, pdfs and cdfs are precomputed for the
            "expected" distance range during connector construction, and looked
            up, rather than recomputed here.

        """
        # Numerator
        if d<self.expected_seq_length:
            numerator = self.stored_pdfs[d]
        else:
            numerator = norm_pdf(d, self._mu, self._sigma)
        
        # Normalize by AUC within the range of observable d values
        max_d = s_dna_len - 1  # Maximum d observable
        # min_d = -1 * max_d  # Minimum d observable
        
        if self._sigma == 0:
            auc = 1.0  # all the gaussian is within the (-(L-1), +(L-1)) range
        else:
            if max_d<self.expected_seq_length:
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
        # The probabity to observe d depends on the number of ways of getting d,
        # given two randomly selected positions within a sequence of lentgth L.
        # There are L - abs(d) ways of getting two random position at
        # distance d within a sequence of length L.
        # The total number of equally likely couples is L*(L-1).
        # Therefore  p(d|null_model) = (L-abs(d)) / (L*(L-1))
        denominator = (s_dna_len - abs(d)) / (s_dna_len * (s_dna_len-1))
        
        # compute additive connector energy term as log-likelihood ratio
        e_connector = np.log2(numerator / denominator)
        
        return e_connector

    def print(self) -> None:
        """Prints the connector mu and sigma values
        """
        print(" m: {} s: {}".format(self._mu, self._sigma))

    def export(self, export_file) -> None:
        """Exports Connector data to the given file

        Args:
            export_file (file): File to export the conector

        """
        export_file.write("\n m: {} s: {}".format(self._mu, self._sigma))

    # pylint: disable=R0201
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
