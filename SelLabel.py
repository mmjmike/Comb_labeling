import tkinter
import csv
import time
import math
import copy

# Globals

RES_TYPES = ("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
             "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
LABEL_TYPES = ("X", "N", "C", "D")
CHECK_PRICE = True
HNCA = False
TIME = False
SEQUENCE = ""
STOCK_TABLE = [[False for col in range(20)]
                    for row in range(4)]
USAGE_TABLE = {}
for res in RES_TYPES:
    USAGE_TABLE.update({res: 1})
USAGE_TABLE["R"] = 2.1
USAGE_TABLE["D"] = 1.8
USAGE_TABLE["E"] = 1.8
USAGE_TABLE["H"] = 1.4
USAGE_TABLE["K"] = 1.2
USAGE_TABLE["M"] = 1.8
USAGE_TABLE["R"] = 1.8
print(USAGE_TABLE)
STOCK_NAME = ""
JOB_NAME = ""
PRICES_TABLE = [[0 for col in range(20)] for row in range(4)]

# Classes


class Sequence:
    """ Sequence class

    it contains the protein sequence and calculates statistics like
    occurrence of each unique pair of residues, ranking of residues by number
    of unique pairs they have and creating the list of residues to be labeled
    with 'calculate stats' method
    based on amino acid stock, written in Stock class
    """

    def __init__(self, name, sequence = ''):
        self.name = name
        self.sequence = sequence
        self.got = True
        if self.sequence == '':
            self.got == False

    def read_from_file(self, filename):

        ## ADD check if file exists

        f = open(filename, 'r')
        seq = f.read()
        f.close()

        ## ADD check input
        ## take expressions from brackets

        self.sequence = seq.upper()
        self.got = True

    def write_to_file(self, filename):
        ## ADD check if file exists

        f = open(filename, 'w')
        f.write(self.sequence)
        f.close()

    def calculate_stats(self, stock):
        # Stock class object required for calculating stats
        # Almost all methods should be rewritten as some of class variables
        # are defined and created in the methods other than '__init__'

        self.stock = stock
        self._rank_residues()
        self._residues_to_label()
        self._count_pairs()
        self._count_nitrogens()
        self._find_bad_residues()
        self._calculate_subtable_coordinates()

    def count_residues(self):
        residue_count = {}
        for residue in RES_TYPES:
            residue_count.update({residue : 0})
        for i in range(len(self.sequence)):
            residue_count[self.sequence[i]] += 1
        return residue_count

    def _rank_residues(self):
        # Bad code for initializing outside __init__

        residue_rank = [0 for _ in RES_TYPES]
        unique_pairs = []
        self.res_no_first = list(RES_TYPES)
        self.res_no_second = list(RES_TYPES)
        for i in range(len(self.sequence) - 1):
            pair = self.sequence[i:i+2]
            if pair[0] in RES_TYPES and pair[1] in RES_TYPES:
                if pair not in unique_pairs:
                    unique_pairs.append(pair)
                    if pair[1] != pair[0]:
                        residue_rank[RES_TYPES.index(pair[1])] += 1
                if pair[0] in self.res_no_first:
                    self.res_no_first.pop(self.res_no_first.index(pair[0]))
                if pair[1] in self.res_no_second:
                    self.res_no_second.pop(self.res_no_second.index(pair[1]))
        self.ranked_residues = [res for rank, res in sorted(zip(residue_rank,
                                                                RES_TYPES))]
        self.ranked_residues.reverse()
        self.rank_of_residue = [rank for rank in sorted(residue_rank)]
        self.rank_of_residue.reverse()
        self.res_no_second.append("P")
        self.residues_first = [res for res in self.ranked_residues if res not in self.res_no_first]
        self.residues_second = [res for res in self.ranked_residues if res not in self.res_no_second]

    def _residues_to_label(self):
        self.residues_nitro = []
        self.residues_not_nitro = []
        self.residues_carbon = []
        self.residues_not_carbon = []
        self.residues_to_label = []
        self.non_labeled_residues = []
        for i in range(len(self.ranked_residues)):
            residue = self.ranked_residues[i]
            if self.rank_of_residue:
                if ('N' in self.stock.label_options[residue]
                        or 'D' in self.stock.label_options[residue]):
                    self.residues_nitro.append(residue)
                else:
                    self.residues_not_nitro.append(residue)
                if ('C' in self.stock.label_options[residue]
                        or 'D' in self.stock.label_options[residue]):
                    self.residues_carbon.append(residue)
                else:
                    self.residues_not_carbon.append(residue)
                if residue in self.residues_carbon or residue in self.residues_nitro:
                    self.residues_to_label.append(residue)
                else:
                    self.non_labeled_residues.append(residue)
        if len(self.residues_not_nitro):
            self.residues_nitro.append("Other")
        if len(self.residues_not_carbon):
            self.residues_carbon.append("Other")

    def _count_nitrogens(self):
        self.min_nitrogens = []
        for i in range(len(self.residues_nitro)):
            pairs_count = 0
            for j in range(len(self.residues_carbon)):
                if self.residue_pairs[j][i]:
                    pairs_count += 1
            self.min_nitrogens.append(pairs_count)

    def _find_bad_residues(self):
        self.bad_residue = []
        for i in range(len(self.residues_nitro)):
            if self.residues_nitro[i] == 'Other':
                self.bad_residue.append(True)
            elif ("D" not in self.stock.label_options[self.residues_nitro[i]]
                  and self.residues_carbon[-1] == 'Other'
                  and self.residues_nitro[i] in self.residues_carbon
                  and self.residue_pairs[self.residues_carbon.index(self.residues_nitro[i])][i]
                  and self.residue_pairs[-1][i]):
                self.bad_residue.append(True)
            else:
                self.bad_residue.append(False)

    def _count_pairs(self):
        self.residue_pairs = [[0 for col in range(len(self.residues_nitro))]
                              for row in range(len(self.residues_carbon))]
        self.all_residue_pairs = [[0 for col in range(len(self.residues_second))]
                                  for row in range(len(self.residues_first))]
        for i in range(len(self.sequence) - 1):  # Count all pairs in sequence,
            pair = self.sequence[i:(i + 2)]
            if pair[0] == '*' or pair[1] == '*':
                continue
            if pair[0] not in self.residues_carbon:
                first_res = len(self.residues_carbon) - 1
            else:
                first_res = self.residues_carbon.index(pair[0])
            if pair[1] not in self.residues_nitro:
                second_res = len(self.residues_nitro) - 1
            else:
                second_res = self.residues_nitro.index(pair[1])
            self.residue_pairs[first_res][second_res] += 1
            index_first = self.residues_first.index(pair[0])
            if pair[1] != "P":
                index_second = self.residues_second.index(pair[1])
                self.all_residue_pairs[index_first][index_second] += 1

    def _calculate_subtable_coordinates(self):
        self.subtable_coordinates = []
        check_other = (self.residues_carbon[-1] == 'Other')
        carbon = 0
        nitro = 0

        for i in range(len(self.residues_to_label)):
            self.new_coordinates = []
            residue = self.residues_to_label[i]
            if i == 0:
                if residue in self.residues_nitro:
                    if residue in self.residues_carbon:
                        if self.residue_pairs[0][0]:
                            self.new_coordinates.append((0, 0))
                        carbon += 1
                    if check_other and self.residue_pairs[-1][0] and not self.bad_residue[0]:
                        self.new_coordinates.append((len(self.residues_carbon) - 1, 0))
                    nitro += 1
            else:
                if residue in self.residues_carbon:
                    for j in range(nitro):
                        if self.residue_pairs[carbon][j]:
                            self.new_coordinates.append((carbon, j))
                    carbon += 1
                if residue in self.residues_nitro:
                    for j in range(carbon):
                        if self.residue_pairs[j][nitro]:
                            self.new_coordinates.append((j, nitro))
                    if check_other and self.residue_pairs[-1][nitro] and not self.bad_residue[nitro]:
                        self.new_coordinates.append((len(self.residues_carbon) - 1, nitro))
                    nitro += 1
            self.subtable_coordinates.append(self.new_coordinates)
        print("Subtable")
        print(self.subtable_coordinates)

    def __str__(self):
        return self.sequence


class Solution:
    '''
    Solution class

    Requires Sequence and Stock objects for initialization

    Has 2 main methods:
        1) Increment state: for each state this method modifies
            self.solution variable depending on whether all the checks
            are met
        2) Check: each solution is checked in each state

    All other methods are supportive for those two
    '''

    def __init__(self, name, sequence, stock):
        self.name = name
        self.sequence = sequence
        self.solution = []
        self.samples_num = 0
        self.found = False
        self.good = False
        self.stock = stock
        self.symmetry = []
        self.depth = 0
        self.price = 0
        self.nitro = 0
        self.carbon = 0
        self.check_pattern_and_pairs = 0
        self.compare_to_previous_code_patterns = 0
        self.check_symmetry = 0
        self.check_subtable = 0
        self.time_increment = 0
        self.time_subtable = 0
        self.time_check = 0
        self.time_nitrogens = 0
        self.time_symmetry = 0
        self.time_patterns = 0
        self.current_patterns = []
        #self.subtable = Subtable()
        self.subtable = []
        self.check_other = (self.sequence.residues_carbon[-1] == 'Other')
        self.labeling_dictionary = {}
        self.best_price = 0
        self.last_price = []
        self.usage = stock.usage
        print(self.usage)

    def __str__(self):
        output = "Solution " + self.name +"\n"
        output += "Patt&pairs\tPrev.patt\tSymmetry\tSubtable\n"
        output += "\t".join(map(str, [self.check_pattern_and_pairs, self.check_symmetry, self.compare_to_previous_code_patterns, self.check_subtable]))
        output += "\n"
        if TIME:
            output += ("time of increment = " + str(self.time_increment)
                       + ". Total increments: " + str(self.check_pattern_and_pairs)
                       + ". Avg: " + str(self.time_increment/self.check_pattern_and_pairs) + "\n")
            output += ("Nitrogens check time = " + str(self.time_nitrogens)
                       + ". Total checks: " + str(self.check_pattern_and_pairs)
                       + ". Avg: " + str(self.time_nitrogens / self.check_pattern_and_pairs) + "\n")
            output += ("Symmetry check time = " + str(self.time_symmetry)
                       + ". Total checks: " + str(self.check_symmetry)
                       + ". Avg: " + str(self.time_symmetry / self.check_symmetry) + "\n")
            output += ("Patterns check time = " + str(self.time_patterns)
                       + ". Total checks: " + str(self.compare_to_previous_code_patterns)
                       + ". Avg: " + str(self.time_patterns / self.compare_to_previous_code_patterns)
                       + "\n")
            output += ("Subtable check time = " + str(self.time_subtable)
                       + ". Total checks: " + str(self.check_subtable)
                       + ". Avg: " + str(self.time_subtable / self.check_subtable)
                       + "\n"
                       )
            output += ("All check time = " + str(self.time_check)
                       + ". Total checks: " + str(self.check_pattern_and_pairs)
                       + ". Avg: " + str(self.time_check / self.check_pattern_and_pairs) + "\n")

        for i in range(len(self.solution)):
            output += self.sequence.residues_to_label[i] + ": " + self.solution[i] + "\n"
        for res in self.sequence.non_labeled_residues:
            output += res + ": " + self._generate_other_code() + "\n"

        if CHECK_PRICE:
            price = self.calculate_total_price()

            output += "\nPrice: " + str(price) + "\n"


        return output

    def calculate_total_price(self):
        price = 0
        for i in range(self.depth):
            residue = self.sequence.residues_to_label[i]
            labeling = self.solution[i]
            price += self.calculate_price_for_residue(residue, labeling)
        return price

    def increment_state(self):
        #t_start = time.time()
        if self.solution == []:

            self.samples_num += 1
            print(self.samples_num)
            self._generate_last_solution()
            self.other_code = self._generate_other_code()
            self.labeling_dictionary.update({'Other': self.other_code})
            self.solution.append(self._base_case())
            self.depth += 1
            if self.sequence.residues_to_label[self.depth - 1] in self.sequence.residues_nitro:
                self.nitro += 1
            if self.sequence.residues_to_label[self.depth - 1] in self.sequence.residues_carbon:
                self.carbon += 1
            self.symmetry.append([self.samples_num])

        elif self.good and self.depth != len(self.sequence.residues_to_label):
            self.symmetry.append(self._calculate_symmetry())
            last_residue = self.sequence.residues_to_label[self.depth - 1]
            last_label = self.solution[-1]
            self.labeling_dictionary.update({last_residue: last_label})
            #self.subtable.add_new_codes(self.new_codes)
            self.subtable.append(self.new_codes)

            if CHECK_PRICE:
                self.last_price.append(self.calculate_price_for_residue(last_residue,
                                                               last_label))
                self.price += self.last_price[-1]

            self.depth += 1

            self.current_patterns.append(self._convert_code_into_nitrogen_pattern(self.solution[-1]))

            #self._extend_subtable()

            self.solution.append(self._base_case())
            if self.sequence.residues_to_label[self.depth - 1] in self.sequence.residues_nitro:
                self.nitro += 1
            if self.sequence.residues_to_label[self.depth - 1] in self.sequence.residues_carbon:
                self.carbon += 1

        elif self.solution[-1] == self.last_solution[self.depth - 1]:

            #print(self)
            if self.sequence.residues_to_label[self.depth - 1] in self.sequence.residues_nitro:
                self.nitro -= 1
            if self.sequence.residues_to_label[self.depth - 1] in self.sequence.residues_carbon:
                self.carbon -= 1
            self.depth -= 1
            self.solution.pop()
            self.symmetry.pop()
            if CHECK_PRICE and self.last_price != []:
                self.price -= self.last_price[-1]
                self.last_price.pop()
            #self.subtable.delete_last_codes()
            if self.subtable != []:
                self.subtable.pop()
            if self.current_patterns != []:
                self.current_patterns.pop()
            self.increment_state()
        else:
            self.solution[-1] = self._increment_code(
                self.solution[-1],
                self.stock.label_options[self.sequence.residues_to_label[self.depth - 1]])


        #self.time_increment += time.time() - t_start

    def _increment_code(self, current_code, options_to_label):
        if len(current_code) == 0:
            return current_code
        if current_code[-1] == options_to_label[-1]:
            return ''.join([self._increment_code(current_code[:-1], options_to_label), options_to_label[0]])
        else:
            return ''.join([current_code[:-1],
                           options_to_label[options_to_label.index(current_code[-1]) + 1]])

    def _base_case(self):
        residue = self.sequence.residues_to_label[len(self.solution)]
        return ''.join([self.stock.label_options[residue][0] for _ in range(self.samples_num)])

    def _calculate_symmetry(self):
        symmetry = self.symmetry[-1]
        code = self.solution[-1]
        if len(symmetry) == self.samples_num:
            return symmetry
        new_symmetry = []
        start_point = 0
        for i in range(len(symmetry)):
            if symmetry[i] > 1:
                number = 1
                for j in range(symmetry[i] - 1):
                    if code[start_point + j] == code[start_point + j + 1]:
                        number += 1
                    else:
                        new_symmetry.append(number)
                        number = 1
                    if j == symmetry[i] - 2:
                        new_symmetry.append(number)
            else:
                new_symmetry.append(1)
            start_point += symmetry[i]
        return new_symmetry

    def check(self):
        #t2_start = time.time()


        result = (self._check_pattern_and_pairs()
                  and self._compare_to_previous_code_patterns()
                  and self._check_symmetry()
                  and self._check_subtable_2()
                  )

        #if len(self.solution) > 3 and self.solution[3] == "CNNNN":
        #    print ("stop")
        #    print(self._check_subtable())
        #    print(self._check_pattern_and_pairs())
        #    print(self._compare_to_previous_code_patterns())


        #if self.depth == len(self.sequence.residues_to_label):
        #
        if result and CHECK_PRICE and self.found:
            result = result and self._check_price()
        self.good = result
        #print(result)

        #self.time_check += time.time() - t2_start
        return result

    def _check_pattern_and_pairs(self):
        #t3_start = time.time()
        #self.check_pattern_and_pairs += 1
        if self.sequence.residues_to_label[self.depth - 1] not in self.sequence.residues_nitro:
            #self.time_nitrogens += time.time() - t3_start
            return True
        nitrogen_factor = 2
        doubles_factor = 2
        if HNCA:
            nitrogen_factor = 3
        nitrogens, doubles = self._count_nitrogens_and_doubles()
        max_pairs = 0
        if nitrogens:
            max_pairs += nitrogen_factor ** nitrogens
        if doubles:
            max_pairs += doubles_factor ** doubles
        return max_pairs >= self.sequence.min_nitrogens[self.nitro - 1]

    def _count_nitrogens_and_doubles(self):
        nitrogens = 0
        doubles = 0
        code = self.solution[-1]
        for i in range(len(code)):
            if code[i] == "N":
                nitrogens += 1
            if code[i] == "D":
                doubles += 1
        return nitrogens, doubles

    def _compare_to_previous_code_patterns(self):
        #t4_start = time.time()
        #self.compare_to_previous_code_patterns += 1
        this_pattern = self._convert_code_into_nitrogen_pattern(self.solution[-1])
        for pattern in self.current_patterns:
            if this_pattern == pattern:
                #self.time_patterns += time.time() - t4_start
                return False
        #self.time_patterns += time.time() - t4_start
        return True

    def _compare_nitrogen_patterns(self, code1, code2):
        return (self._convert_code_into_nitrogen_pattern(code1)
                == self._convert_code_into_nitrogen_pattern(code2))

    def _convert_code_into_nitrogen_pattern(self, code):
        new_code = ''
        for i in range(len(code)):
            if code[i] == 'C' or code[i] == 'X':
                new_code += '0'
            if code[i] == 'N' or code[i] == 'D':
                new_code += '1'
        return new_code

    def _check_subtable(self):
        #t1_start = time.time()

        #self.check_subtable += 1

        self.new_codes = []
        current_residue = self.sequence.residues_to_label[self.depth - 1]
        residues_carbon = self.sequence.residues_carbon
        residues_nitro = self.sequence.residues_nitro
        current_labeling = self.solution[-1]
        residue_pairs = self.sequence.residue_pairs
        self.labeling_dictionary.update({self.sequence.residues_to_label[self.depth - 1]: self.solution[-1]})

        if self.depth == 1:
            if current_residue in residues_nitro:
                if current_residue in residues_carbon and residue_pairs[0][0]:
                    code_1 = self._calculate_code_multiple(current_labeling,
                                                                   current_labeling)
                    self.new_codes.append(code_1)
                if (self.check_other and residue_pairs[-1][0] and not self.sequence.bad_residue[0]):
                    code_2 = self._calculate_code_multiple(self._generate_other_code(),
                                                                   current_labeling)
                    if residue_pairs[0][0] and code_2 == code_1:
                        #self.time_subtable += time.time() - t1_start
                        return False
                    self.new_codes.append(code_2)
        else:
            if current_residue in residues_carbon:
                carbon_index = residues_carbon.index(current_residue) ## self.carbon - 1
                for i in range(self.nitro):
                    if residue_pairs[carbon_index][i] and residues_nitro[i] != current_residue:
                        second_residue_code = self.labeling_dictionary[residues_nitro[i]]
                        new_code = self._calculate_code_multiple(current_labeling,
                                                                       second_residue_code)
                        if self._code_used(new_code):
                            #self.time_subtable += time.time() - t1_start
                            return False
                        self.new_codes.append(new_code)
            if current_residue in residues_nitro:
                nitro_index = residues_nitro.index(current_residue)   ## self.nitro - 1
                for i in range(self.carbon):
                    if residue_pairs[i][nitro_index]:
                        first_residue_code = self.labeling_dictionary[residues_carbon[i]]
                        new_code = self._calculate_code_multiple(first_residue_code,
                                                                       current_labeling)
                        if self._code_used(new_code):
                            #self.time_subtable += time.time() - t1_start
                            return False
                        self.new_codes.append(new_code)
                if (self.check_other and residue_pairs[-1][nitro_index] and not self.sequence.bad_residue[nitro_index]):
                    new_code = self._calculate_code_multiple(self._generate_other_code(),
                                                                   current_labeling)
                    if self._code_used(new_code):
                        #self.time_subtable += time.time() - t1_start
                        return False
                    self.new_codes.append(new_code)
        #self.time_subtable += time.time() - t1_start
        return True

    def _check_subtable_2(self):
        #t1_start = time.time()

        #self.check_subtable += 1

        self.new_codes = set()
        self.labeling_dictionary.update({self.sequence.residues_to_label[self.depth - 1]: self.solution[-1]})
        coordinates_list = self.sequence.subtable_coordinates[self.depth - 1]

        for coordinates in coordinates_list:
            first_residue = self.sequence.residues_carbon[coordinates[0]]
            second_residue = self.sequence.residues_nitro[coordinates[1]]
            first_residue_labeling = self.labeling_dictionary[first_residue]
            second_residue_labeling = self.labeling_dictionary[second_residue]
            code = self._calculate_code_multiple(first_residue_labeling,
                                                 second_residue_labeling)
            if self._code_used(code):
                return False
            self.new_codes.add(code)
        return True

    def _code_used(self, new_code):
        if new_code in self.new_codes:
            return True
        for code_list in self.subtable:
            if new_code in code_list:
                return True
        return False

    def _check_symmetry(self):
        #t5_start = time.time()
        #self.check_symmetry += 1
        symmetry = self.symmetry[-1]
        code = self.solution[-1]
        options = self.stock.label_options[self.sequence.residues_to_label[self.depth - 1]]
        if len(symmetry) == self.samples_num:
            #self.time_symmetry += time.time() - t5_start
            return True
        start_point = 0
        for i in range(len(symmetry)):
            if symmetry[i] > 1:
                for j in range(symmetry[i] - 1):
                    if (options.index(code[start_point + j])
                            > options.index(code[start_point + j + 1])):
                        #self.time_symmetry += time.time() - t5_start
                        return False
            start_point += symmetry[i]
        #self.time_symmetry += time.time() - t5_start
        return True

    def _check_price(self):
        last_residue = self.sequence.residues_to_label[self.depth - 1]
        last_label = self.solution[-1]

        price = (self.price
                 + self.calculate_price_for_residue(last_residue, last_label)
                 + self._predict_price()
                 )
        if price >= self.best_price:
            return False
        return True

    def _predict_price(self):
        price = 0
        dict_filled = {}
        for i in range(len(self.sequence.residues_to_label) - self.depth):
            residue = self.sequence.residues_to_label[i + self.depth]
            dict_filled.update({residue: 0})
        for i in range(len(self.sequence.min_nitrogens) - self.nitro):
            # add comparison between N and D variants
            residue = self.sequence.residues_nitro[i + self.nitro]
            if residue == "Other":
                continue
            num_of_nitrogens = math.ceil(math.log(self.sequence.min_nitrogens[i + self.nitro])
                               / math.log(2))
            dict_filled[residue] += num_of_nitrogens
            if 'N' in self.stock.label_options[residue]:
                price += (num_of_nitrogens
                          * int(self.stock.price_dict[residue]["N"])
                          * self.usage[residue])
            else:
                price += (num_of_nitrogens
                          * int(self.stock.price_dict[residue]["D"])
                          * self.usage[residue])
        for i in range(len(self.sequence.residues_carbon) - self.carbon):
            # add comparison between N and D variants
            residue = self.sequence.residues_carbon[i + self.carbon]
            if residue == "Other":
                continue
            num_of_carbons = 0

            for j in range(len(self.sequence.residues_nitro)):
                if self.sequence.residue_pairs[i + self.carbon][j] and self.sequence.residue_pairs[-1][j] and self.check_other:
                    num_of_carbons = 1
                    break
            dict_filled[residue] += num_of_carbons
            if 'C' in self.stock.label_options[residue]:
                price += (num_of_carbons
                          * int(self.stock.price_dict[residue]["C"])
                          * self.usage[residue])
            else:
                price += (num_of_carbons
                          * int(self.stock.price_dict[residue]["D"])
                          * self.usage[residue])
                if 'N' not in self.stock.label_options[residue]:
                    dict_filled[residue] += num_of_carbons
        for i in range(len(self.sequence.residues_to_label) - self.depth):
            residue = self.sequence.residues_to_label[i + self.depth]
            cheapest_option = self._find_cheapest_option(residue)
            price += (cheapest_option
                      * (self.samples_num - dict_filled[residue]) \
                      * self.usage[residue])
        return price

    def _find_cheapest_option(self, residue):
        first = True
        cheapest = 0
        for option in self.stock.label_options[residue]:
            if first:
                cheapest = int(self.stock.price_dict[residue][option])
                first = False
            elif int(self.stock.price_dict[residue][option]) < cheapest:
                cheapest = int(self.stock.price_dict[residue][option])
        return cheapest

    def _calculate_code_multiple(self, first_aa_code, second_aa_code):
        samples = self.samples_num
        code_string = ''.join([str(self._calculate_code(first_aa_code[i],
                                                        second_aa_code[i])) for i in range(samples)])
        return code_string

    def _calculate_code(self, first_aa_code, second_aa_code):
        if HNCA:
            if second_aa_code == 'N':
                if first_aa_code == 'X' or first_aa_code == 'N':
                    code = 1
                elif first_aa_code == 'C':
                    code = 2
                else:
                    code = 3
            elif second_aa_code == 'D':
                if first_aa_code == 'X' or first_aa_code == 'N':
                    code = 4
                else:
                    code = 3
            else:
                code = 0
        else:
            if second_aa_code == 'N' or second_aa_code == 'D':
                if first_aa_code == 'C' or first_aa_code == 'D':
                    code = 2
                else:
                    code = 1
            else:
                code = 0
        return code

    def calculate_price_for_residue(self, residue, labeling):
        price = 0
        for i in range(len(labeling)):
            price += int(self.stock.price_dict[residue][labeling[i]]) * self.usage[residue]
        return price

    def _generate_last_solution(self):
        self.last_solution = []
        for residue in self.sequence.residues_to_label:
            current_code = ''.join([self.stock.label_options[residue][-1] for _ in range(self.samples_num)])
            self.last_solution.append(current_code)

    def _generate_other_code(self):
        code = ''
        for i in range(self.samples_num):
            code += 'X'
        return code

    def print_table(self):
        aa_types_1 = self.sequence.residues_carbon
        aa_types_2 = self.sequence.residues_nitro
        aa_pairs = self.sequence.residue_pairs
        code_dict = self.labeling_dictionary
        num_of_samples = self.samples_num

        cell_height = 20
        cell_width = 50



        row_number = len(aa_types_1)
        column_number = len(aa_types_2)

        top = tkinter.Tk()
        C = tkinter.Canvas(top, bg="black", height=cell_height * (row_number + 2),
                           width=cell_width * (column_number + 2))


        ### Draw the table and fill the headers
        for i in range(row_number):
            C.create_line(0, cell_height * (i + 1), cell_width * (column_number + 2), cell_height * (i + 1),
                          fill="white")
            C.create_text(cell_width / 2, cell_height * (i + 2.5), text=aa_types_1[i], fill="white")
            C.create_text(cell_width * 1.5, cell_height * (i + 2.5), text=code_dict[aa_types_1[i]], fill="white")
        for i in range(column_number):
            C.create_line(cell_width * (i + 1), 0, cell_width * (i + 1), cell_height * (row_number + 2),
                          fill="white")
            C.create_text(cell_width * (i + 2.5), cell_height / 2, text=aa_types_2[i], fill="white")

            column_codes = self.get_codes(aa_types_2, code_dict)
            codes_match = 0
            for k in range(column_number):
                if (k != i and self._compare_nitrogen_patterns(column_codes[i], column_codes[k])):
                    codes_match = 1
                if codes_match:
                    cell_color = 'red'
                else:
                    cell_color = 'green'
            C.create_rectangle(cell_width * (i + 2) + 1, cell_height + 1, cell_width * (i + 3) - 1,
                               cell_height * 2 - 1, fill=cell_color)

            C.create_text(cell_width * (i + 2.5), cell_height * 1.5, text=code_dict[aa_types_2[i]], fill="white")
        C.create_line(0, cell_height * (row_number + 1), cell_width * (column_number + 2),
                      cell_height * (row_number + 1), fill="white")
        C.create_line(cell_width * (column_number + 1), 0, cell_width * (column_number + 1),
                      cell_height * (row_number + 2), fill="white")

        all_aa_codes = [['' for col in range(column_number)] for row in range(row_number)]
        for i in range(row_number):
            for j in range(column_number):
                if aa_pairs[i][j]:
                    all_aa_codes[i][j] = self._calculate_code_multiple(code_dict[aa_types_1[i]], code_dict[aa_types_2[j]])

        ### fill the cells
        good_cells = 0
        bad_cells = 0
        for i in range(row_number):
            for j in range(column_number):
                if aa_pairs[i][j] > 0:
                    codes_match = 0
                    for k in range(row_number):
                        if (k != i and all_aa_codes[i][j] == all_aa_codes[k][j]):
                            codes_match = 1
                    if codes_match:
                        cell_color = 'red'
                        bad_cells += 1
                    else:
                        cell_color = 'green'
                        good_cells += 1
                    C.create_rectangle(cell_width * (j + 2) + 1, cell_height * (i + 2) + 1,
                                       cell_width * (j + 3) - 1, cell_height * (i + 3) - 1, fill=cell_color)
                    # C.create_text(cell_width * (j + 2.5), cell_height * (i + 2.5), text = aa_pairs[i][j], fill="white")
                    C.create_text(cell_width * (j + 2.5), cell_height * (i + 2.5), text=all_aa_codes[i][j],
                                  fill="white")
        C.create_rectangle(0, 0, cell_width * 2 - 1, cell_height * 2 - 1, fill="black")
        print(str(good_cells / (good_cells + bad_cells) * 100) + '%', (good_cells + bad_cells))

        C.pack()
        top.mainloop()

    def get_codes(self, aa_string, code_dict):
        aa_codes = []
        for i in range(len(aa_string)):
            aa_codes.append(code_dict[aa_string[i]])
        return aa_codes


class Stock:
    '''
Stock object

Contains the amino acid stock
and also prices for each type of labeling
'''

    def __init__(self, name):
        self.name = name
        self.label_dict = {}
        self.label_options = {} ##Подумать, может свойством сделать?
        self.price_dict = {}
        self.usage = {}
        for res in RES_TYPES:
            self.usage.update({res: 1})

    def read_from_file(self, filename):

        # ADD check if file exists

        with open(filename, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            d = list(reader)

        # ADD check if file format is correct

        for i in range(len(d) - 1):
            self.label_dict.update({d[i + 1][0]: ''.join(d[i + 1][1:])})
        self._generate_label_options()

    def read_from_table(self, table):
        for i in range(len(RES_TYPES)):
            res = RES_TYPES[i]
            table_part = [table[j][i] for j in range(len(table))]
            self.label_dict.update({res: ''.join(["1" if cell else "0" for cell in table_part])})
        self._generate_label_options()

    def read_label_prices(self, input_file):
        with open(input_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            d = list(reader)

        for i in range(len(d) - 1):
            residue = d[i+1][0]
            residue_prices = {
                'X': d[i+1][1],
                'N': d[i+1][2],
                'C': d[i+1][3],
                'D': d[i+1][4]
            }
            self.price_dict.update({residue: residue_prices})

    def read_prices_from_table(self, table):
        pass

    def add_usage_dict(self, dict):
        self.usage = dict

    def _generate_label_options(self):
        for residue in RES_TYPES:
            if residue in self.label_dict:
                stock_to_change = self.label_dict[residue]
                option = []
                ## Labeling types are in preferred order,
                ##but due to success with this software the order doesn't matter.
                if stock_to_change[1] == '1':
                    option.append('N')
                if stock_to_change[2] == '1':
                    option.append('C')
                if stock_to_change[3] == '1':
                    option.append('D')
                if stock_to_change[0] == '1':
                    option.append('X')
                self.label_options.update({residue : option})

    def read_prices(self, filename):

        # ADD check if file exists

        with open(filename, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            d = list(reader)
            label_types = d[0][1:]

        # ADD check if file format is correct

        for i in range(len(d) - 1):
            curr_dict = {}
            for j in range(label_types):
                curr_dict.update({label_types[j] : int(d[i + 1][j + 1])})
            self.price_dict.update({d[i + 1][0] : curr_dict})

    def write_to_file(self, filename):

        #ADD check if file exists

        f = open(filename, 'w')
        f.write('Res\t' + '\t'.join(LABEL_TYPES) + '\n')
        for residue in RES_TYPES:
            line = residue
            for label in LABEL_TYPES:
                line += '\t'
                if residue in self.price_dict:
                    line += str(self.price_dict[residue][label])
                else:
                    line += '0'
            line += '\n'
            f.write(line)
        f.close()

    def print(self):
        i = 1
        for key in self.label_dict:
            print (str(i) + "." + key + ":" + self.label_dict[key])
            i += 1
        print("\n")

    def write_prices(self, filename):
        ## RAW METHOD!!!!!!!!!!!!!!!!
        # ADD check if file exists

        f = open(filename, 'w')
        f.write('Res\t' + '\t'.join(LABEL_TYPES) + '\n')
        for residue in RES_TYPES:
            line = residue
            for label in LABEL_TYPES:
                line += '\t'
                if (residue in self.label_options
                    and label in self.label_options[residue]):
                    line += '1'
                else:
                    line += '0'
            line += '\n'
            f.write(line)
        f.close()


class Subtable:

    def __init__(self):
        self.table_of_tables = []
        self.map = set()

    def add_new_codes(self, list_of_codes):
        self.table_of_tables.append(list_of_codes)
        self._recalculate_map()

    def delete_last_codes(self):
        if self.table_of_tables != []:
            self.table_of_tables.pop()
            self._recalculate_map()

    def map_code(self, code):
        return code in self.map

    def _recalculate_map(self):
        self.map = set()
        for table in self.table_of_tables:
            self.map = self.map.union(table)


class SolutionFinder:

    def __init__(self, options):
        self.name = options["job_name"]
        self.sequence = Sequence(self.name, options["sequence"])
        self.stock = Stock(self.name, options["label_table"])
        self.stock.add_prices(options["prices_table"])
        self.sequence.calculate_stats(self.stock)
        self.solution = Solution(self.name, self.sequence, self.stock)
        self.check_price = options["optimize_price"]
        self.hnca = options["HNCA"]
        self.paused = False

    def run(self):
        pass

    def pause(self):
        self.paused = True

    def resume(self):
        self.paused = False

    def stop(self):
        pass

    def check_state(self):
        pass

# FUNCTIONS

def set_parameters(parameters):
    global SEQUENCE, HNCA, CHECK_PRICE, STOCK_TABLE, STOCK_NAME, JOB_NAME, PRICES_TABLE
    JOB_NAME = parameters["job_name"]
    SEQUENCE = parameters["sequence"]
    STOCK_TABLE = parameters["label_table"]
    HNCA = parameters["HNCA"]
    PRICES_TABLE = parameters["prices_table"]
    CHECK_PRICE = parameters["optimize_price"]

def main():

    # DEFINITION OF OBJECTS: Stock, Sequence, Solution

    # stock = Stock(JOB_NAME)
    # stock.read_from_table(STOCK_TABLE)
    # if CHECK_PRICE:
    #     stock.read_prices_from_table(PRICES_TABLE)
    # sequence = Sequence(JOB_NAME, sequence=SEQUENCE)
    # sequence.calculate_stats(stock)
    # solution = Solution(JOB_NAME, sequence, stock)

    stock = Stock("Lyukmanova")
    stock.read_from_file("LabelStore.txt")
    stock.read_label_prices("prices.txt")

    stock2 = Stock("Goncharuk")
    stock2.read_from_file("LabelStore_3_ot_seregi.txt")
    stock2.read_label_prices("prices.txt")
    stock2.add_usage_dict(USAGE_TABLE)

    stock3 = Stock("Kesha")
    stock3.read_from_file("LabelStore_Maslennikov.txt")
    stock3.read_label_prices("prices.txt")

    sequence = Sequence("Kv2.1")
    sequence.read_from_file("kv21.seq")
    sequence.calculate_stats(stock)

    sequence2 = Sequence("KdpD")
    sequence2.read_from_file("KdpD.seq")
    sequence2.calculate_stats(stock3)

    sequence3 = Sequence("Nav1D")
    sequence3.read_from_file("nav1D.seq")
    sequence3.calculate_stats(stock2)

    solution = Solution("Nav1D", sequence3, stock2)


    # END OF DEFINITION

    # Variables

    samples_number = 0 # if we check price - this var is to remember at what samples number
                       # solution was found, so we stop searching the next number of samples
    iteration = 0 # just count iterations (optional)
    t0 = time.time()
    solutions = 0
    final_depth = len(solution.sequence.residues_to_label)
    ten_minute_counter = 1
    best_solution = copy.deepcopy(solution)

    while True:

        #   MAIN CYCLE, increments Solution and checks it

        iteration += 1
        solution.increment_state()
        solution.check()

        # Check if solution is found
        if solution.depth == final_depth and solution.good:

            # block that works if we optimize the solution by price
            if CHECK_PRICE:

                if not solution.found:
                    samples_number = solution.samples_num  # remember samples number
                    solution.best_price = solution.calculate_total_price()
                if solution.best_price > solution.calculate_total_price():
                    best_solution = copy.deepcopy(solution)
                    solution.best_price = solution.calculate_total_price()
                    print("curr best price: " + str(solution.best_price))
                    print(solution)
                solution.found = True
            else:
                solution.found = True
                best_solution = copy.deepcopy(solution)
                break
            solutions += 1

            #print("iteration: " + str(iteration))
            #print("solutions: " + str(solutions))
            #print(solution)

        # break the cycle if number of samples increased and the solution is already found
        # bad idea to continue the search anyways
        if solution.found and samples_number != solution.samples_num:
            break

        # TESTING BLOCK

        if iteration % 1000000 == 0:
            print("iteration: " + str(iteration) + "; time elapsed: " + str(round(time.time() - t0, 2)))
            print("iterations/min:", str(round(iteration * 60 / (time.time() - t0), 1)))
            if CHECK_PRICE:
                print("curr best price: " + str(solution.best_price))
            print("")
            #if not solution.found:
                #print(solution)

        if time.time() - t0 > 600 * ten_minute_counter:
            print(str(ten_minute_counter * 10) + " minutes passed")
            print("iteration: " + str(iteration) + "; time elapsed: " + str(round(time.time() - t0, 2)))
            if CHECK_PRICE:
                print("curr best price: " + str(solution.best_price))
            print("")
            ten_minute_counter += 1

        #if time.time() - t0 > 1800:
         #   break

        ## END OF TESTING BLOCK
    print(iteration)
    print("Best solution:")
    print(best_solution)
    print("solutions: " + str(solutions))

    t1 = time.time()
    print(t1 - t0)

    return best_solution


if __name__ == '__main__':
    solution = main()
    solution.print_table()