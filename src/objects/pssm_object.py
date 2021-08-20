"""P object
   PSSM object is a type of recognizer implementing sequence-specific recognition.
   The PSSM object stores all PSSM data and its configuration.
   
   A PSSM (Position-Specific Scoring Matrix) is a model of the sequence-specific
   recognition of DNA by a protein. Given a DNA sequence of length L, a PSSM of
   length L returns a scalar representing the binding energy provided by the
   recognition. This binding energy is derived from the log-likelihood between
   the model's preference for the bases in the given sequence (at their given
   location) and the null hypothesis (here, for generality, that all bases are
   equally probable).
"""

import random
import numpy as np
import decimal as dec


class PssmObject():
    """PSSM object is a type of recognizer object
    """

    def __init__(self, pwm, config: dict) -> None:
        """PSSM constructor:
            - Gets/sets:
                - the length of the PSSM
                - the frequency matrix
                - all the PSSM-specific config parameters
            

        Args:
            pwm (numpy.array): PWM
            config: configuration from JSON file
        """
        
        # set PSSM length and position weight matrix
        self.length = len(pwm)  # length of the numpy array
        self.pwm = pwm  # numpy array of dictionaries
        self.pssm = None #scoring matrix
        
        # assign PSSM-specific configuration elements
        self.mutate_probability_random_col = config[
            "MUTATE_PROBABILITY_RANDOM_COL"
        ]
        self.mutate_probability_mutate_col = config[
            "MUTATE_PROBABILITY_MUTATE_COL"
        ]
        self.mutate_probability_flip_cols = config[
            "MUTATE_PROBABILITY_FLIP_COL"
        ]
        self.mutate_probability_flip_rows = config[
            "MUTATE_PROBABILITY_FLIP_ROW"
        ]
        self.mutate_probability_shift_left = config[
            "MUTATE_PROBABILITY_SHIFT_LEFT"
        ]
        self.mutate_probability_shift_right = config[
            "MUTATE_PROBABILITY_SHIFT_RIGHT"
        ]
        self.mutate_probability_increase_pwm = config[
            "MUTATE_PROBABILITY_INCREASE_PWM"
                ]
        self.mutate_probability_decrease_pwm = config[
            "MUTATE_PROBABILITY_DECREASE_PWM"
                ]
        self.min_columns = config["MIN_COLUMNS"]
        self.max_columns = config["MAX_COLUMNS"]

        self.pseudo_count = config["PSEUDO_COUNT"]

        self.upper_print_probability = config["UPPER_PRINT_PROBABILITY"]
        self.scan_reverse_complement = config["SCAN_REVERSE_COMPLEMENT"]
        
        # Compute PSSM Matrix based on PWM
        self.recalculate_pssm()


    def update_length(self):
        """Updates the length attribute.
        """
        self.length = len(self.pwm)


    def mutate(self, org_factory) -> None:
        """Mutation operators associated to the PSSM recognizer

        Args:
            org_factory (OrganismFactory): Creates objects
        """
        
        # Code to keep track if mutations that shift the boundaries of the
        # PSSM placement have occurred (shift left/shift right/increase pwm/
        # decrease pwm).
        pssm_displacement_code = [0, 0]
        
        if random.random() < self.mutate_probability_random_col:

            # Randomize PSSM column 
            # [substitutes column with a randomized column]
            new_col = org_factory.get_pwm_column()
            # Select a random col in self.pwm
            column_to_update = random.randint(0, self.length - 1)
            # Insert it in that position
            self.pwm[column_to_update] = new_col
        
        if random.random() < self.mutate_probability_mutate_col:
            '''Mutate a column of the PSSM. The mutation involves transferring
            some of the weight from a "donor" base to an "acceptor" base'''
            # Chose a random column of that PSSM
            idx_of_random_col = random.randint(0, self.length - 1)
            
            # Chose from which base to which base the transfer of weight will occur
            donor_base, acceptor_base = random.sample(['a','c','g','t'], 2)
            
            # Current values of the two bases
            donor_current_prob = self.pwm[idx_of_random_col][donor_base]            
            acceptor_current_prob = self.pwm[idx_of_random_col][acceptor_base]
            
            no_BSs = org_factory.pwm_number_of_binding_sites
            donor_current_count = donor_current_prob * no_BSs            
            acceptor_current_count = acceptor_current_prob * no_BSs
            
            # entity of the transfer of counts
            transfer = random.randint(0, int(donor_current_count) - 1)
            
            # New count values for the two bases
            donor_new_count = donor_current_count - transfer
            acceptor_new_count = acceptor_current_count + transfer
            
            # Back to probabilities
            donor_new_prob = donor_new_count / no_BSs
            acceptor_new_prob = acceptor_new_count / no_BSs
            
            # To avoid extra decimals due to floating point errors:
            # Define the number of decimal digits actually required.
            # All frequencies are multiples of 1/no_BSs
            smallest_freq = dec.Decimal('1') / dec.Decimal(str(no_BSs))
            # Number of decimals in the smallest frequency
            no_decimals = len(str(smallest_freq).split(".")[1])
            # No frequency value needs more decimal digits than  smallest_freq.
            # Therefore we can round according to  no_decimals
            donor_new_prob = round(donor_new_prob, no_decimals)
            acceptor_new_prob = round(acceptor_new_prob, no_decimals)
            # Update pwm
            self.pwm[idx_of_random_col][donor_base] = donor_new_prob
            self.pwm[idx_of_random_col][acceptor_base] = acceptor_new_prob

        if random.random() < self.mutate_probability_flip_cols:
            # Swaps two PSSM columns
            # col1 --> col2, col2 --> col1
            col1, col2 = random.sample(range(self.length), 2)
            # Select two random columns and swap them
            tmp_col = self.pwm[col1]
            self.pwm[col1] = self.pwm[col2]
            self.pwm[col2] = tmp_col

        if random.random() < self.mutate_probability_flip_rows:
            # Swaps two PSSM rows
            
            # Pick rows (A, C, T or G) to be swapped
            bases = ["a", "c", "g", "t"]
            random.shuffle(bases)
            base1, base2 = bases[:2]

            # Swap rows
            for i in range(self.length):
                tmp_base = self.pwm[i][base1]
                self.pwm[i][base1] = self.pwm[i][base2]
                self.pwm[i][base2] = tmp_base

        if random.random() < self.mutate_probability_shift_left:
            # Shift PSSM from right to left, rolling over
            self.pwm = np.roll(self.pwm, 1)
            
            # The left bound of the PSSM shifts 1 bp to the left (-1)
            pssm_displacement_code[0] -= 1
            # The right bound of the PSSM shifts 1 bp to the left (-1)
            pssm_displacement_code[1] -= 1

        if random.random() < self.mutate_probability_shift_right:
            # Shift PSSM from left to right, rolling over
            self.pwm = np.roll(self.pwm, -1)
            
            # The left bound of the PSSM shifts 1 bp to the right (+1)
            pssm_displacement_code[0] += 1
            # The right bound of the PSSM shifts 1 bp to the right (+1)
            pssm_displacement_code[1] += 1

        if random.random() < self.mutate_probability_increase_pwm:
            # Increase length of PWM
            if self.length < self.max_columns: # do only if allowed

                #generate a new column
                new_col = org_factory.get_pwm_column()

                # Add the new column to one side (chose randomly left or right)
                if random.random() < 0.5:
                    # Insert to the left
                    tmp_array = [new_col] + self.pwm.tolist()
                    
                    # The left bound of the PSSM shifts 1 bp to the left (-1)
                    pssm_displacement_code[0] -= 1
                    
                else:
                    # Insert to the right
                    tmp_array = self.pwm.tolist() + [new_col]
                    
                    # The right bound of the PSSM shifts 1 bp to the right (+1)
                    pssm_displacement_code[1] += 1

                # assign newly made PWM
                self.pwm = np.array(tmp_array)
                # Update length attribute
                self.update_length()

        if random.random() < self.mutate_probability_decrease_pwm:
            # Decrease length of PWM
            if self.length > self.min_columns:
                
                # Remove a column from one side (chose randomly left or right)
                if random.random() < 0.5:
                    # Remove from the left
                    self.pwm = self.pwm[1:]
                    
                    # The left bound of the PSSM shifts 1 bp to the right (+1)
                    pssm_displacement_code[0] += 1
                
                else:
                    # Remove from the right
                    self.pwm = self.pwm[:-1]
                    
                    # The right bound of the PSSM shifts 1 bp to the left (-1)
                    pssm_displacement_code[1] -= 1
                
                # Update length attribute
                self.update_length()
        
        # recompute PSSM
        # mutation operators affect the PWM (frequency matrix)
        # so the PSSM is re-computed after mutations take place
        self.recalculate_pssm()
        
        # If the PSSM boundaries have changed, report the pssm-displacement
        # code, so that connectors can eventually be adjusted if necessary
        if pssm_displacement_code != [0, 0]:
            return pssm_displacement_code

    # Calculate self.pssm based on self.pwm
    def recalculate_pssm(self) -> None:
        """ Calculates the PSSM based on the pwm values
        """
        tmp_pssm = []
        for column in self.pwm:
            # From pwm to pssm
            # log2(base/0.25) = log2(4.0*base)
            decimals = 2
            tmp_bases = []
            # cast to float so round function does not become crazy
            tmp_bases.append(
                float(np.log2(4.0 * column["c"] + self.pseudo_count))
            )
            tmp_bases.append(
                float(np.log2(4.0 * column["t"] + self.pseudo_count))
            )
            tmp_bases.append(
                float(np.log2(4.0 * column["g"] + self.pseudo_count))
            )
            tmp_bases.append(
                float(np.log2(4.0 * column["a"] + self.pseudo_count))
            )

            tmp_pssm.append(
                {
                    "c": round(tmp_bases[0], decimals),
                    "t": round(tmp_bases[1], decimals),
                    "g": round(tmp_bases[2], decimals),
                    "a": round(tmp_bases[3], decimals),
                }
            )
        # Assign re-computed PSSM
        self.pssm = np.array(tmp_pssm)


    def get_score(self, s_dna: str) -> float:
        """Get the score for given DNA secuence (s_dna)

        Args:
            s_dna: dna partial sequence (length of the pssm)

        Returns:
            score assigned to s_dna. If reverse sequence is better, reverse
            score is returned
        """

        complement = {"a": "t", "t": "a", "g": "c", "c": "g"}
        # gets a score from pssm
        score = 0
        score_reverse = float("-inf")
        str_length = len(s_dna)
      
        # score given strand
        for i in range(str_length):
            score += self.pssm[i][s_dna[i]]

        # if reverse sequence scoring is activated, score reverse
        if self.scan_reverse_complement:
            score_reverse = 0
            for i in range(str_length):
                score_reverse += self.pssm[str_length - i - 1][
                                 complement[s_dna[str_length - i - 1]]]
               
        # Return the max binding score
        # (only truly applies if scan_reverse_complement is activated)
        if score_reverse > score:
            return(score_reverse)
        else:
            return(score)
    

    def print(self) -> None:
        """Print PSSM object (similar to Logo format)
           Prints a consensus sequence, with uppercase characters
           depending on user-defined threshold: upper_print_probability
        """

        recognized = ""

        for position in self.pwm:
            base = "a"
            # Find max base
            for new_base in position.keys():
                if position[new_base] > position[base]:
                    base = new_base

            if position[base] >= self.upper_print_probability:
                base = base.upper()
            recognized += base

        print(recognized)

    def export(self, export_file) -> None:
        """Exports pssm to a file

        Args:
            export_file: File to write the output
        """
        recognized = ""

        for position in self.pwm:
            base = "a"
            # Find max base
            for new_base in position.keys():
                if position[new_base] > position[base]:
                    base = new_base
            # Change to uppercase based on probability
            if position[base] >= self.upper_print_probability:
                base = base.upper()
            recognized += base
        
        export_file.write("\n" + recognized)

    def is_connector(self) -> bool:
        """node is not a connector

        Returns:
            False because is a pssm recognizer
        """
        return False

    def is_pssm(self):
        """node is a pssm recognizer

        Returns:
            True because is a pssm recognizer
        """
        return True
