import tkinter
import csv
import time
import math
import copy
import residues as rs


#
# Вывод таблицы кодов и расшифровки - типы спектров номер кода
# Опция проверки решения
# Предсказание цен: разобраться с азотами и углеродами
# Вырожденные метки, придумать, как с ними быть
# Причесать все по PEP8
# Вбить реальные прайсы
# Мануал
# Подготовить входные файлы
# Написать кучу комментариев

# Globals

RND_DEBUG = True
ONLY_NUMBERS = True
CONFIG_FILE = "comblabel.config"
RES_TYPES = ("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
             "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
LABEL_TYPES = ("N", "C", "D", "X", "A", "S", "T", "F")
CN_TYPES = ("D", "S", "T")
C_TYPES = ("C", "D", "A", "S", "T", "F")

CHECK_PRICE = False
HNCA = False
LOHR = False
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

#print(USAGE_TABLE)
STOCK_NAME = ""
JOB_NAME = ""
PRICES_TABLE = [[0 for col in range(20)] for row in range(4)]

SPEC_VEC = [1, 1, 1, 1, 1, 1]
LABEL_VEC = [1, 1, 1, 1, 1, 0, 0, 0]




# Classes


class Sequence:
    """ Sequence class

    it contains the protein sequence and calculates statistics like
    occurrence of each unique pair of residues, ranking of residues by number
    of unique pairs they have and creating the list of residues to be labeled
    with 'calculate stats' method
    based on amino acid stock, written in Stock class
    """

    def __init__(self, task, name="", sequence=""):
        self.name = name
        self.sequence = sequence
        self.got = True
        self.task = task
        self.nitro_types = ("N", "D", "S", "T")
        self.carbon_types = ("C", "D", "A", "S", "T", "F")

    def read_from_file(self, filename):

        ## ADD check if file exists

        sequence_lines = get_lines(filename)
        if sequence_lines == 1:
            print("ERROR! Sequence file '" + filename + "' not found")
            return
        new_lines = []
        if [0][0] == ">":
            for line in sequence_lines[1:]:
                if line[0] == ">":
                    break
                new_lines.append(line)
        else:
            new_lines = sequence_lines
        sequence = clean_sequence(new_lines)
        if len(sequence) < 2:
            print("ERROR! Sequence too short or absent")
            return
        ## ADD check input
        ## take expressions from brackets

        self.sequence = sequence.upper()
        self.got = True
        return 0

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
        self._check_carbon_pairs()
        self._count_nitrogens()
        self._find_bad_residues()
        self._calculate_subtable_coordinates()
        print(self.residues_to_label)

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
        # self.ranked_residues = [res for rank, res in sorted(zip(residue_rank,
        #                                                         RES_TYPES))]
        # self.ranked_residues.reverse()
        # self.rank_of_residue = [rank for rank in sorted(residue_rank)]
        # self.rank_of_residue.reverse()
        self.ranked_residues, self.rank_of_residue = self._sort_residues(RES_TYPES, residue_rank)
        # print(self.rank_of_residue)
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
            if self.rank_of_residue[i]:
                nitro = False
                for type in self.nitro_types:
                    if type in self.stock.label_options[residue]:
                        nitro = True
                        break
                if nitro:
                    self.residues_nitro.append(residue)
                else:
                    self.residues_not_nitro.append(residue)

                carbon = False
                for type in self.carbon_types:
                    if type in self.stock.label_options[residue]:
                        carbon = True
                        break
                if carbon:
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

    def _sort_residues(self, residue_list, residue_rank):
        residues = list(residue_list)
        print(residues)
        for i in range(len(residues)-1):
            for j in range(len(residues)-1-i):
                if residue_rank[i] < residue_rank[i+j+1]:
                    temp_res = residues[i]
                    temp_rank = residue_rank[i]
                    residue_rank[i] = residue_rank[i + j + 1]
                    residues[i] = residues[i + j + 1]
                    residue_rank[i + j + 1] = temp_rank
                    residues[i + j + 1] = temp_res
        return residues, residue_rank

    def _check_carbon_pairs(self):
        carb_other = 0
        changed = False

        # print('here', self.residues_carbon)
        if self.residues_carbon and self.residues_carbon[-1] == "Other":
            carb_other = 1
        if self.residues_nitro and self.residues_nitro[-1] == "Other":
            carbons = len(self.residues_carbon) - carb_other
            for i in range(carbons):
                k = carbons - 1 - i
                got_pair = False
                for j in range(len(self.residues_nitro) - 1):
                    if self.residue_pairs[k][j]:
                        got_pair = True
                        break
                if not got_pair:
                    self.residues_not_carbon.append(self.residues_carbon[k])
                    residue = self.residues_carbon[k]
                    del self.residues_carbon[k]
                    if residue not in self.residues_nitro:
                        del self.residues_to_label[self.residues_to_label.index(residue)]
                    changed = True
        self.ranks = {}
        for res in self.residues_to_label:
            self.ranks.update({res: 0})
        if changed:
            self._count_pairs()
        for i in range(len(self.residues_nitro)):
            if self.residues_nitro[i] == "Other":
                break
            for j in range(len(self.residues_carbon)):
                if self.residue_pairs[j][i]:
                    if self.residues_carbon[j] != "Other":
                        self.ranks[self.residues_carbon[j]] += 1
                    if self.residues_carbon[j] != self.residues_nitro[i]:
                        self.ranks[self.residues_nitro[i]] += 1
        if RND_DEBUG:
            ranks = []
            residues = []
            for residue in self.residues_to_label:
                if residue in self.ranks:
                    residues.append(residue)
                    ranks.append(self.ranks[residue])
            self.residues_to_label, self.residue_ranks = self._sort_residues(residues, ranks)
        else:
            self.residues_to_label = list(sorted(self.ranks, key=self.ranks.get))
            self.residues_to_label.reverse()
        print(self.residues_to_label)
        self.residues_nitro_new = []
        self.residues_carbon_new = []
        for res in self.residues_to_label:
            if res in self.residues_carbon:
                self.residues_carbon_new.append(res)
            if res in self.residues_nitro:
                self.residues_nitro_new.append(res)
        if "Other" in self.residues_carbon:
            self.residues_carbon_new.append("Other")
        if "Other" in self.residues_nitro:
            self.residues_nitro_new.append("Other")
        self.residues_carbon = self.residues_carbon_new
        self.residues_nitro = self.residues_nitro_new
        self._count_pairs()

    def _add_rank(self, residue):
        self.ranks[self.residues_to_label.index(residue)] += 1

    def _count_nitrogens(self):
        self.min_nitrogens = []
        for i in range(len(self.residues_nitro)):
            pairs_count = 0
            for j in range(len(self.residues_carbon)):
                if self.residue_pairs[j][i]:
                    pairs_count += 1
            self.min_nitrogens.append(pairs_count)

    def _find_bad_residues(self):
        self.bad_residues = []
        for i in range(len(self.residues_nitro)):
            res = self.residues_nitro[i]
            if res == "Other":
                continue
            if ("D" not in self.stock.label_options[res]
                  and "S" not in self.stock.label_options[res]
                  and "T" not in self.stock.label_options[res]
                  and self.residues_carbon[-1] == 'Other'
                  and self.residues_nitro[i] in self.residues_carbon
                  and self.residue_pairs[self.residues_carbon.index(res)][i]
                  and self.residue_pairs[-1][i]):
                self.bad_residues.append(res)

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

        # carbon_total = len(self.residues_carbon)
        # nitro_total = len(self.residues_nitro)
        # if self.residues_nitro[-1] == "Other":
        #     nitro_total -= 1
        # total_cells = carbon_total * nitro_total
        total_cells = 0

        for i in range(len(self.residues_to_label)):
            self.new_coordinates = []
            residue = self.residues_to_label[i]
            if i == 0:
                if residue in self.residues_nitro:
                    if residue in self.residues_carbon:
                        if self.residue_pairs[0][0]:
                            self.new_coordinates.append((0, 0))
                        carbon += 1
                    if check_other and self.residue_pairs[-1][0] and residue not in self.bad_residues:
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
                    if check_other and self.residue_pairs[-1][nitro] and residue not in self.bad_residues:
                        self.new_coordinates.append((len(self.residues_carbon) - 1, nitro))
                    nitro += 1
            total_cells += len(self.new_coordinates)
            self.subtable_coordinates.append(self.new_coordinates)
        curr_cells = 0
        print("Total cells: ", total_cells)
        for i in range(len(self.subtable_coordinates)):
            add = len(self.subtable_coordinates[i])
            curr_cells += add
            print("Curr cells: ", curr_cells,". Added: ", add)
            # if curr_cells > total_cells / 2:
            #     print("Median at: ", i+1)
            #     break

    def _need_carbon(self, residue):
        carbon_number = self.residues_carbon.index(residue)
        if self.residues_carbon[-1] == "Other":
            scan_nitro = len(self.residues_nitro)
            if self.residues_nitro[-1] == "Other":
                scan_nitro -= 1
            if residue in self.residues_nitro:
                nitro_index = self.residues_nitro.index(residue)
                if self.residue_pairs[carbon_number][nitro_index] and self.residue_pairs[-1][nitro_index] and residue not in self.bad_residues:
                    return self._find_cheapest_option(residue, CN_TYPES)
            for i in range(scan_nitro):
                if self.residue_pairs[carbon_number][i] and self.residue_pairs[-1][i]:
                    if residue == self.residues_nitro[i]:
                        if residue in self.bad_residues:
                            continue
                        else:
                            return self._find_cheapest_option(residue, CN_TYPES)
                    else:
                        return self._find_cheapest_option(residue, self.carbon_types)
        return self._find_cheapest_option(residue, RES_TYPES)

    def _find_cheapest_option(self, residue, labeling_options):
        options = self.stock.label_options[residue]
        first_option = True
        lower_price = 0
        cheapest_option = ""
        for i in range(len(options)):
            option = options[i]
            if option in labeling_options:
                curr_price = self.stock.price_dict[residue][option]
                if first_option:
                    lower_price = curr_price
                    cheapest_option = option
                    first_option = False
                elif curr_price < lower_price:
                    cheapest_option = option
                    lower_price = curr_price
                else:
                    continue
        return cheapest_option

    def _need_nitrogen(self, residue, res_options):
        min_nitrogens = self.min_nitrogens[self.residues_nitro.index(residue)]
        all_options = self.stock.label_options[residue]
        label_powers = self.task.coding_table.label_power
        nitro_options = []
        self.variants = []
        current_power = 1
        for option in all_options:
            if option in self.nitro_types and label_powers[option] > 1:
                nitro_options.append(option)
        if nitro_options == []:
            return 1
        for option in res_options:
            if option in self.nitro_types:
                current_power *= label_powers[option]
        if (min_nitrogens == 1 and current_power > 1) or (min_nitrogens > 1 and current_power >= min_nitrogens):
            return []
        else:
            self._add_nitrogen([], nitro_options, residue, min_nitrogens, current_power)
            best_result = []
            first_iteration = True
            for variant in self.variants:
                if first_iteration:
                    first_iteration = False
                    best_result = variant
                    best_price = self._calc_price_for_variant(residue, variant)
                else:
                    if self._calc_price_for_variant(residue, variant) < best_price:
                        best_result = variant
                        best_price = self._calc_price_for_variant(residue, variant)
            return best_result

    def _calc_price_for_variant(self, residue, variant):
        price = 0
        for type in variant:
            price += self.stock.price_dict[residue][variant]
        return price

    def _add_nitrogen(self, list1, nitro_options, residue, min_nitrogens, current_power):
        for option in nitro_options:
            list1.append(option)
            if self._check_power(list1, min_nitrogens, current_power):
                self.variants.append(list1)
            else:
                self._add_nitrogen(list1, nitro_options, residue, min_nitrogens, current_power)
            list1.pop()

    def _check_power(self, list1, min_nitrogens, current_power):
        power = current_power
        for option in list1:
            power *= self.task.coding_table.label_power[option]
        return power > min_nitrogens

    def _calculate_base_price(self, residue, options):
        price = 0
        for option in options:
            price += self.stock.price_dict[residue][option]
        return price

    def predict_prices(self):
        residues = self.residues_to_label
        self.required_options = []
        self.base_prices = []
        self.add_prices = []
        self.max_requirement = 0
        for res in residues:
            res_options = []
            if res in self.residues_nitro:
                if res in self.residues_carbon:
                    res_options.append(self._need_carbon(res))
                nitro_required = self._need_nitrogen(res, res_options)
                if nitro_required == 1:
                    return 1
                for nitro in nitro_required:
                    res_options.append(nitro)
            else:
                res_options.append(self._need_carbon(res))
            self.required_options.append(res_options)
            if len(res_options) > self.max_requirement:
                self.max_requirement = len(res_options)
        for i in range(len(residues)):
            res = residues[i]
            cheapest_option = self._find_cheapest_option(res, RES_TYPES)
            add_price = self.stock.price_dict[res][cheapest_option]
            base_price = self._calculate_base_price(res, self.required_options[i])
            for j in range(self.max_requirement - len(self.required_options[i])):
                base_price += add_price
            self.base_prices.append(base_price)
            self.add_prices.append(add_price)
        return 0

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

    def __init__(self, task):
        self.task = task
        self.name = self.task.name
        self.sequence = self.task.sequence
        self.coding_table = self.task.coding_table
        self.optimize_price = self.task.optimize_price
        if self.optimize_price:
            self.base_prices = self.sequence.base_prices
            self.add_prices = self.sequence.add_prices
            self.max_requirement = self.sequence.max_requirement
            self.predicted_prices = [0 for res in self.sequence.residues_to_label]
        self.solution = []
        self.samples_num = 0
        self.found = False
        self.good = False
        self.stock = self.task.stock
        self.symmetry = []
        self.depth = 0
        self.max_depth = 0
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
        self.usage = self.stock.usage
        self.rt = rs.make_rt()
        self.codes_dict = self.coding_table.codes_dict
        self.label_power = self.coding_table.label_power

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
        output += ("Max depth: " + str(self.max_depth))
        if self.optimize_price:
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

    def read_existing_solution(self, filename):
        # rewrite, so that it can read self-made solutions
        self.labeling_dictionary = {}
        with open(filename, 'r') as f:
            reader = csv.reader(f)
            d = list(reader)

            # ADD check if file format is correct
        for i in range(len(d)):
            self.labeling_dictionary.update({d[i][0]: d[i][1]})
        self.samples_num = len(d[0][1])
        self.labeling_dictionary.update({"Other": self._generate_other_code()})
        print(self.labeling_dictionary)

    def increment_state(self):
        #t_start = time.time()
        if self.solution == []:

            self.samples_num += 1
            print("Current max depth: ", self.max_depth)
            print(self.samples_num)
            if self.optimize_price:
                if self.samples_num >= self.max_requirement:
                    self.predicted_prices = [(self.samples_num - self.max_requirement)*x+y for x,y in zip(self.add_prices, self.base_prices)]
                else:
                    self.predicted_prices = self.base_prices
            self._generate_last_solution()
            self.other_code = self._generate_other_code()
            self.labeling_dictionary.update({'Other': self.other_code})
            self.solution.append(self._base_case())
            self.depth += 1
            if self.depth > self.max_depth:
                self.max_depth = self.depth
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

            if self.optimize_price:
                self.last_price.append(self.calculate_price_for_residue(last_residue,
                                                               last_label))
                self.price += self.last_price[-1]

            self.depth += 1
            if self.depth > self.max_depth:
                self.max_depth = self.depth
            if last_residue in self.sequence.residues_nitro:
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
                if self.current_patterns != []:
                    self.current_patterns.pop()
            if self.sequence.residues_to_label[self.depth - 1] in self.sequence.residues_carbon:
                self.carbon -= 1
            self.depth -= 1
            self.solution.pop()
            self.symmetry.pop()
            if self.optimize_price and self.last_price != []:
                self.price -= self.last_price[-1]
                self.last_price.pop()
            #self.subtable.delete_last_codes()
            if self.subtable != []:
                self.subtable.pop()

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
                  #and self._compare_to_previous_code_patterns()
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
        if result and self.optimize_price and self.found:
            result = result and self._check_price()
        self.good = result
        #print(result)


        #self.time_check += time.time() - t2_start
        return result

    def _check_pattern_and_pairs(self):
        #t3_start = time.time()
        #self.check_pattern_and_pairs += 1
        current_residue = self.sequence.residues_to_label[self.depth - 1]
        if current_residue not in self.sequence.residues_nitro:
            #self.time_nitrogens += time.time() - t3_start
            return True

        ## Переписать для более общего случая:
        else:
            max_pairs = 1
            code = self.solution[-1]
            for labeling in code:
                max_pairs *= self.label_power[labeling]
        # nitrogen_factor = self.label_power["N"]
        # doubles_factor = self.label_power["D"]
        # nitrogens, doubles = self._count_nitrogens_and_doubles()
        # max_pairs = 0
        # if nitrogens:
        #     max_pairs += (nitrogen_factor ** nitrogens)
        # if doubles:
        #     max_pairs *= (doubles_factor ** doubles)
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
        if self.sequence.residues_to_label[self.depth - 1] not in self.sequence.residues_nitro:
            return True
        if LOHR:
            this_pattern = self._convert_code_into_ncd_pattern(self.solution[-1])
        else:
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

    def _compare_ncd_patterns(self, code1, code2):
        return (self._convert_code_into_ncd_pattern(code1)
                == self._convert_code_into_ncd_pattern(code2))

    def _convert_code_into_nitrogen_pattern(self, code):
        new_code = ''
        for i in range(len(code)):
            if code[i] == 'C' or code[i] == 'X':
                new_code += '0'
            if code[i] == 'N' or code[i] == 'D':
                new_code += '1'
        return new_code

    def _convert_code_into_ncd_pattern(self, code):
        new_code = ''
        for i in range(len(code)):
            if code[i] == 'C' or code[i] == 'X' or code[i] == 'A':
                new_code += '0'
            if code[i] == 'N':
                new_code += '1'
            if code[i] == 'D':
                new_code += '2'
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
                if (self.check_other and residue_pairs[-1][0] and current_residue not in self.bad_residues):
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
                if (self.check_other and residue_pairs[-1][nitro_index] and current_residue not in self.bad_residues):
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
        last_residue = self.sequence.residues_to_label[self.depth-1]
        last_label = self.solution[-1]
        predicted_price = 0
        for i in range(len(self.solution) - self.depth):
            predicted_price += self.predicted_prices[self.depth+i]

        price = (self.price
                 + self.calculate_price_for_residue(last_residue, last_label)
                 + predicted_price
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
                          * rs.rt.one_let[residue].usage)
            else:
                price += (num_of_nitrogens
                          * int(self.stock.price_dict[residue]["D"])
                          * rs.rt.one_let[residue].usage)
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
                          * rs.rt.one_let[residue].usage)
            else:
                price += (num_of_carbons
                          * int(self.stock.price_dict[residue]["D"])
                          * rs.rt.one_let[residue].usage)
                if 'N' not in self.stock.label_options[residue]:
                    dict_filled[residue] += num_of_carbons
        for i in range(len(self.sequence.residues_to_label) - self.depth):
            residue = self.sequence.residues_to_label[i + self.depth]
            cheapest_option = self._find_cheapest_option(residue)
            price += (cheapest_option
                      * (self.samples_num - dict_filled[residue])
                      * rs.rt.one_let[residue].usage)
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
        return self.codes_dict[first_aa_code][second_aa_code]

    def calculate_price_for_residue(self, residue, labeling):
        price = 0
        for i in range(len(labeling)):
            price += self.stock.price_dict[residue][labeling[i]]
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

    def print_table(self, codes=True):
        aa_types_1 = self.sequence.residues_carbon
        aa_types_2 = self.sequence.residues_nitro
        aa_pairs = self.sequence.residue_pairs
        code_dict = self.labeling_dictionary
        num_of_samples = self.samples_num
        other_nitro = 0
        if aa_types_2[-1] == "Other":
            other_nitro = 1

        cell_height = 20
        cell_width = 50



        row_number = len(aa_types_1)
        column_number = len(aa_types_2) - other_nitro

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
                    if codes:
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
                    else:
                        cell_color = 'green'
                        good_cells += 1
                    C.create_rectangle(cell_width * (j + 2) + 1, cell_height * (i + 2) + 1,
                                       cell_width * (j + 3) - 1, cell_height * (i + 3) - 1, fill=cell_color)
                    if codes:
                        C.create_text(cell_width * (j + 2.5), cell_height * (i + 2.5), text=all_aa_codes[i][j],
                                      fill="white")
                    else:
                        C.create_text(cell_width * (j + 2.5), cell_height * (i + 2.5), text = aa_pairs[i][j], fill="white")

        C.create_rectangle(0, 0, cell_width * 2 - 1, cell_height * 2 - 1, fill="black")
        print(str(good_cells / (good_cells + bad_cells) * 100) + '%', (good_cells + bad_cells))

        C.pack()
        top.mainloop()

    def get_codes(self, aa_string, code_dict):
        aa_codes = []
        for i in range(len(aa_string)):
            aa_codes.append(code_dict[aa_string[i]])
        return aa_codes

    def make_codes_table(self, filename):
        self.codes_table = []
        self.meanings_table = []
        for res in self.sequence.non_labeled_residues:
            self.labeling_dictionary[res] = self._generate_other_code()
        for i in range(len(self.sequence.sequence) - 1):
            pair = self.sequence.sequence[i:i+2]
            residue1 = pair[0]
            residue2 = pair[1]
            scheme1 = self.labeling_dictionary[residue1]
            scheme2 = self.labeling_dictionary[residue2]
            code = self._calculate_code_multiple(scheme1, scheme2)
            meaning = "".join([residue1, str(i+1), " - ", residue2, str(i+2)])
            self.codes_table.append(code)
            self.meanings_table.append(meaning)

        self.ranked_meanings = [mean for code, mean in sorted(zip(self.codes_table,
                                                                self.meanings_table))]
        self.ranked_codes_table = [code for code in sorted(self.codes_table)]
        with open(filename, 'w') as f:
            for i in range(len(self.ranked_codes_table)):
                f.write(": ".join([self.ranked_codes_table[i], self.ranked_meanings[i]]) + "\n")
            f.close()

    def write_to_file(self, filename):
        output = "Res"
        for i in range(self.samples_num):
            output += ",sample_" + str(i+1)
        output += "\n"
        for i in range(len(self.solution)):
            output += self.sequence.residues_to_label[i] + "," + ",".join(list(self.solution[i])) + "\n"
        for res in self.sequence.non_labeled_residues:
            output += res + "," + ",".join(list(self._generate_other_code())) + "\n"

        with open(filename, 'w') as f:
            f.write(output)
            f.close()

    def write_full_pairs_table(self, filename):
        output = "," + ",".join(self.sequence.residues_second) + "\n"
        for i in range(len(self.sequence.residues_first)):
            output += self.sequence.residues_first[i]
            for j in range(len(self.sequence.residues_second)):
                # MAYBE REPLACE 0 BY SOME SYMBOL????
                output += "," + str(self.sequence.all_residue_pairs[i][j])
            if i+1 < len(self.sequence.residues_first):
                output += "\n"
        with open(filename, 'w') as f:
            f.write(output)
            f.close()

    def write_pairs_table(self, filename):
        output = ","
        if self.sequence.residues_nitro[-1] == "Other":
            output += ",".join(self.sequence.residues_nitro[:-1])
            output += "," + "".join(self.sequence.residues_not_nitro)
        else:
            output += ",".join(self.sequence.residues_nitro)
        output += "\n"
        for i in range(len(self.sequence.residues_carbon)):
            res1 = self.sequence.residues_carbon[i]
            if res1 == "Other":
                output += "".join(self.sequence.residues_not_carbon)
            else:
                output += res1
            for j in range(len(self.sequence.residues_nitro)):
                # MAYBE REPLACE 0 BY SOME SYMBOL????
                output += "," + str(self.sequence.residue_pairs[i][j])
            if i + 1 < len(self.sequence.residues_carbon):
                output += "\n"
        with open(filename, 'w') as f:
            f.write(output)
            f.close()

    def write_pairs_codes(self, filename):
        output = ",,"
        if self.sequence.residues_nitro[-1] == "Other":
            res_list = self.sequence.residues_nitro[:-1]
        else:
            res_list = self.sequence.residues_nitro
        output += ",".join(res_list)
        output += "\n,"
        for res in res_list:
            output += "," + self.labeling_dictionary[res]
        output += "\n"
        for i in range(len(self.sequence.residues_carbon)):
            res1 = self.sequence.residues_carbon[i]
            if res1 == "Other":
                output += "".join(self.sequence.residues_not_carbon)
            else:
                output += res1
            output += "," + self.labeling_dictionary[res1]
            for j in range(len(self.sequence.residues_nitro)):
                res2 = self.sequence.residues_nitro[j]
                if self.sequence.residue_pairs[i][j]:
                    code = self._calculate_code_multiple(self.labeling_dictionary[res1],
                                                         self.labeling_dictionary[res2])
                else:
                    code = ""
                output += "," + code
            if i + 1 < len(self.sequence.residues_carbon):
                output += "\n"
        with open(filename, 'w') as f:
            f.write(output)
            f.close()


class Stock:
    '''
Stock object

Contains the amino acid stock
and also prices for each type of labeling
'''

    def __init__(self, name=""):
        self.name = name
        self.label_dict = {}
        self.label_options = {} ##Подумать, может свойством сделать?
        self.price_dict = {}
        self.usage = {}
        self.label_types = []
        self.price_label_types = []
        for res in RES_TYPES:
            self.usage.update({res: 1})

    def read_from_file(self, filename):

        try:
            with open(filename, 'r') as f:
                reader = csv.reader(f, delimiter='\t')
                d = list(reader)
        except FileNotFoundError:
            print("ERROR! Stock file '" + filename + "' not found")
            return 1



        # ADD check if file format is correct
        try:
            self.label_types = d[0][1:]
            for i in range(len(d) - 1):
                self.label_dict.update({d[i + 1][0]: ''.join(d[i + 1][1:])})
            self._generate_label_options()
        except IndexError:
            print("ERROR in stock file '" + filename + "'.\nPlease check the length of each row.\nUse Tab as separator")
            return 1
        return 0

    def read_from_table(self, table):
        #self.label_types = ['X', 'N', 'C', 'D']
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
            residue_prices = {}
            for j in range(len(self.label_types)):
                residue_prices.update({self.label_types[j]: d[i+1][j+1]})
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
                for label in LABEL_TYPES:
                    if label in self.label_types:
                        label_index = self.label_types.index(label)
                        if stock_to_change[label_index] == '1':
                            option.append(label)
                self.label_options.update({residue: option})

    def read_prices(self, filename):

        # ADD check if file exists
        lines = get_lines(filename)
        if lines == 1:
            print("ERROR! Stock file '" + filename + "' not found")
            return 1
        else:
            d = []
            for line in lines:
                if line[0] != "#":
                    d.append(line.split(","))
            self.price_label_types = d[0][1:]
            for label in self.label_types:
                if label not in self.price_label_types:
                    print("ERROR! Prices are not specified for '" + label + "' label type")
                    return 1
            print(d)

            # ADD check if file format is correct
            try:
                for i in range(len(d) - 1):
                    curr_dict = {}
                    for j in range(len(self.price_label_types)):
                        try:
                            price = float(d[i + 1][j + 1])
                        except ValueError:
                            print("ERROR! Price must be set in digits")
                            return 1
                        curr_dict.update({self.price_label_types[j]: price})
                    residue_type = d[i + 1][0]
                    print("curr dict: ", self.label_options[residue_type])
                    for label_type in self.label_options[residue_type]:
                        if curr_dict[label_type] < 0:
                            print("ERROR! Price is not specified or negative for '"
                                  + residue_type + "'-residue's '"
                                  + label_type + "' label type")
                            return 1
                    self.price_dict.update({residue_type: curr_dict})
                print("Here", self.price_dict)
            except IndexError:
                print("ERROR in stock file '" + filename + "'.\nPlease check the length of each row.\nUse Tab as separator")
                return 1
            return 0

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


class Task:

    def __init__(self):
        self.sequence = Sequence(self)
        self.stock = Stock()
        self.optimize_price = False
        self.spectra_vector = []
        self.labeling_types_vector = []
        self.set = False
        self.name = ""
        self.sequence_ok = False
        self.stock_ok = False
        self.price_ok = False
        self.spectra_ok = False
        self.label_types_ok = False

    def import_sequence(self, sequence):
        self.sequence.sequence = sequence
        self.sequence_ok = True

    def read_stock(self, stock_file):
        result = self.stock.read_from_file(stock_file)
        if result == 0:
            self.stock_ok = True
        return result

    def read_spec_vec(self, spec_vec):
        if len(spec_vec) != 7:
            print("ERROR! Not for all (7) spectra the usage is specified")
            return 1
        for item in spec_vec:
            if item != "1" and item != "0":
                print("ERROR! Only '1' or '0' must be used to specify spectra")
                return 1
        self.spectra_vector = [int(item) for item in spec_vec]
        self.spectra_ok = True
        return 0

    def read_label_vec(self, lab_vec):
        if len(lab_vec) != 8:
            print("ERROR! Not for all (8) labeling types specified")
            return 1
        for item in lab_vec:
            if item != "1" and item != "0":
                print("ERROR! Only '1' or '0' must be used to specify labeling types")
                return 1
        self.labeling_types_vector = [int(item) for item in lab_vec]
        self.label_types_ok = True
        return 0

    def read_sequence(self, seq_file):
        result = self.sequence.read_from_file(seq_file)
        if result == 0:
            self.sequence_ok = True
        return result

    def set_price_optimization(self, value):
        self.optimize_price = value
        self.price_ok = not value

    def read_prices(self, price_file):
        result = self.stock.read_prices(price_file)
        if result == 0:
            self.price_ok = True
        return result

    def create_solution(self):
        self.sequence.calculate_stats(self.stock)
        self.coding_table = CodingTable(self)
        if self.optimize_price:
            if self.sequence.predict_prices():
                return 1
        self.solution = Solution(self)
        return 0

    def solve(self):
        # print("sequence: ", self.sequence.sequence, "\n")
        # print("stock: \n")
        # self.stock.print()
        # print("Spectra: ", self.spectra_vector, "\n")
        # print("Label types: ", self.labeling_types_vector, "\n")
        # print("Optimize price: ", self.optimize_price)
        if self.sequence_ok and self.stock_ok and self.price_ok and self.spectra_ok and self.label_types_ok:
            if self.create_solution():
                return 1
            samples_number = 0
            iteration = 0
            t0 = time.time()
            solutions = 0
            final_depth = len(self.solution.sequence.residues_to_label)
            ten_minute_counter = 1
            self.best_solution = copy.deepcopy(self.solution)

            while True:

                #   MAIN CYCLE, increments Solution and checks it

                iteration += 1
                self.solution.increment_state()
                self.solution.check()

                # Check if solution is found
                if self.solution.depth == final_depth and self.solution.good:

                    # block that works if we optimize the solution by price
                    if self.optimize_price:

                        if not self.solution.found:
                            samples_number = self.solution.samples_num  # remember samples number
                            self.solution.best_price = self.solution.calculate_total_price()
                            self.best_solution = copy.deepcopy(self.solution)
                        if self.solution.best_price > self.solution.calculate_total_price():
                            self.best_solution = copy.deepcopy(self.solution)
                            self.solution.best_price = self.solution.calculate_total_price()
                            print("curr best price: " + str(self.solution.best_price))
                            print(self.solution)
                        self.solution.found = True
                    else:
                        self.solution.found = True
                        self.best_solution = copy.deepcopy(self.solution)
                        break
                    solutions += 1

                    # print("iteration: " + str(iteration))
                    # print("solutions: " + str(solutions))
                    # print(solution)

                # break the cycle if number of samples increased and the solution is already found
                # bad idea to continue the search anyways
                if self.solution.found and samples_number != self.solution.samples_num:
                    break

                # TESTING BLOCK

                if iteration % 1000000 == 0:
                    print("iteration: " + str(iteration) + "; time elapsed: " + str(round(time.time() - t0, 2)))
                    print("iterations/min:", str(round(iteration * 60 / (time.time() - t0), 1)))
                    if self.optimize_price:
                        print("curr best price: " + str(self.solution.best_price))
                    print("")
                    print(self.solution)
                    # if not solution.found:
                    # print(solution)

                if time.time() - t0 > 600 * ten_minute_counter:
                    print(str(ten_minute_counter * 10) + " minutes passed")
                    print("iteration: " + str(iteration) + "; time elapsed: " + str(round(time.time() - t0, 2)))
                    if CHECK_PRICE:
                        print("curr best price: " + str(self.solution.best_price))
                    print("")
                    ten_minute_counter += 1
                    print(self.solution)

                    # if time.time() - t0 > 1800:
                    #   break

                    ## END OF TESTING BLOCK
            # print("Iteration: ", iteration)
            # print("Best solution:")
            # print(self.best_solution)
            # print("solutions: " + str(solutions))

            t1 = time.time()
            # print(self.stock.label_options)
            # print(t1 - t0)
            self.best_solution.print_table(not ONLY_NUMBERS)
            return iteration

    def write_results(self):
        solution_filename = self.name + "_solution.txt"
        code_dict_filename = self.name + "_code_dictionary.txt"
        codes_filename = self.name + "_codes.txt"
        pairs_table_filename = self.name + "_pairs.txt"
        full_table_filename = self.name + "_all_pairs.txt"
        pairs_codes_table_filename = self.name + "_pairs_codes.txt"
        self.best_solution.write_to_file(solution_filename)
        self.best_solution.make_codes_table(code_dict_filename)
        self.coding_table.write_to_file(codes_filename)
        self.best_solution.write_full_pairs_table(full_table_filename)
        self.best_solution.write_pairs_table(pairs_table_filename)
        self.best_solution.write_pairs_codes(pairs_codes_table_filename)


class CodingTable:

    def __init__(self, task):
        self.task = task
        self.label_vec = task.labeling_types_vector
        self.spec_vec = task.spectra_vector
        self.spec_types = ("HSQC", "HNCO", "HNCA", "HNCOCA", "CO-HNCA", "DQ-HNCA", "HNCACO")
        self.spec_list = [self._HSQC, self._HNCO, self._HNCA, self._HNCOCA, self._COHNCA, self._DQHNCA, self._HNCACO]
        self.label_meanings = {
            "X": "000",
            "N": "100",
            "C": "001",
            "D": "111",
            "A": "010",
            "S": "101",
            "T": "110",
            "F": "011",
        }
        self.label_types = ("N", "C", "D", "X", "A", "S", "T", "F")  # THIS order of labels sets their priority in search for solution
        self.letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K']
        self.codes_dict = {}
        self.label_power = {}
        self.vectors = []
        self.spectra_numbers = []
        self._make_coding_table()

    def _HSQC(self, atom_list):
        return atom_list[3]

    def _HNCO(self, atom_list):
        return atom_list[3] and atom_list[2]

    def _HNCA(self, atom_list):
        return atom_list[3] and (atom_list[1] or atom_list[4])

    def _HNCOCA(self, atom_list):
        return atom_list[3] and atom_list[2] and atom_list[1]

    def _COHNCA(self, atom_list):
        return atom_list[3] and atom_list[1] and not atom_list[2]

    def _DQHNCA(self, atom_list):
        return atom_list[3] and atom_list[1] and atom_list[4]

    def _HNCACO(self, atom_list):
        return atom_list[3] and atom_list[4] and atom_list[5]

    def _make_coding_table(self):

        label_types_list = []
        for i in range(len(self.label_types)):
            if self.label_vec[i]:
                label_types_list.append(self.label_types[i])
        spec_types_list = []
        for i in range(len(self.spec_list)):
            if self.spec_vec[i]:
                spec_types_list.append(self.spec_list[i])
                self.spectra_numbers.append(i)
        codes_table = [[0 for i in range(len(label_types_list))] for j in range(len(label_types_list))]
        self.vectors = [[0 for _ in spec_types_list]]
        for i in range(len(label_types_list)):
            label_2 = label_types_list[i]
            for j in range(len(label_types_list)):
                label_1 = label_types_list[j]
                atom_list = self._make_atom_list(label_1, label_2)
                vector = []
                for spectrum in spec_types_list:
                    vector.append(spectrum(atom_list))
                if vector in self.vectors:
                    code = self.vectors.index(vector)
                else:
                    code = len(self.vectors)
                    self.vectors.append(vector)
                if code > 9:
                    codes_table[j][i] = self.letters[code-10]
                else:
                    codes_table[j][i] = str(code)
        for i in range(len(label_types_list)):
            result = []
            for row in codes_table:
                result.append(row[i])
            power = len(set(result))
            self.label_power[label_types_list[i]] = power

        for i in range(len(label_types_list)):
            label_1 = label_types_list[i]
            subdict = {}
            for j in range(len(label_types_list)):
                label_2 = label_types_list[j]
                subdict[label_2] = codes_table[i][j]
            self.codes_dict[label_1] = subdict
        # print(self.codes_dict)

    def _make_atom_list(self, first_type, second_type):
        atom_string = self.label_meanings[first_type] + self.label_meanings[second_type]
        return [int(symbol) for symbol in atom_string]

    def write_to_file(self, filename):
        output = "Code"
        for i in range(len(self.spectra_numbers)):
            output += "," + self.spec_types[self.spectra_numbers[i]]
        output += "\n"
        for i in range(len(self.vectors)):
            if i < 10:
                code = str(i)
            else:
                code = self.letters[i-10]
            output += code + ","
            output += ",".join([str(item) for item in self.vectors[i]])
            if i+1 < len(self.vectors):
                output += "\n"

        with open(filename, 'w') as f:
            f.write(output)
            f.close()


# FUNCTIONS

def set_parameters(parameters):
    global SEQUENCE, HNCA, CHECK_PRICE, STOCK_TABLE, STOCK_NAME, JOB_NAME, PRICES_TABLE
    JOB_NAME = parameters["job_name"]
    SEQUENCE = parameters["sequence"]
    STOCK_TABLE = parameters["label_table"]
    HNCA = parameters["HNCA"]
    PRICES_TABLE = parameters["prices_table"]
    CHECK_PRICE = parameters["optimize_price"]


def main2():

    # DEFINITION OF OBJECTS: Stock, Sequence, Solution

    # stock = Stock(JOB_NAME)
    # stock.read_from_table(STOCK_TABLE)
    # if CHECK_PRICE:
    #     stock.read_prices_from_table(PRICES_TABLE)
    # sequence = Sequence(JOB_NAME, sequence=SEQUENCE)
    # sequence.calculate_stats(stock)
    # solution = Solution(JOB_NAME, sequence, stock)

    # stock = Stock("Lyukmanova")
    # stock.read_from_file("LabelStore.txt")
    # stock.read_label_prices("prices.txt")
    #
    # stock2 = Stock("Goncharuk")
    # stock2.read_from_file("LabelStore_3_ot_seregi.txt")
    # stock2.read_label_prices("prices.txt")
    # stock2.add_usage_dict(USAGE_TABLE)
    #
    # stock3 = Stock("Kesha")
    # stock3.read_from_file("LabelStore_Maslennikov.txt")
    # stock3.read_label_prices("prices.txt")
    # stock3.add_usage_dict(USAGE_TABLE)

    stock4 = Stock("Kesha")
    stock4.read_from_file("LabelStore_Lohr_2.txt")
    # stock4.read_label_prices("prices.txt")
    # stock4.add_usage_dict(USAGE_TABLE)

    # sequence = Sequence("Kv2.1")
    # sequence.read_from_file("kv21.seq")
    # sequence.calculate_stats(stock)
    #
    # sequence2 = Sequence("KdpD")
    # sequence2.read_from_file("KdpD.seq")
    # sequence2.calculate_stats(stock3)
    #
    # sequence3 = Sequence("Nav1D")
    # sequence3.read_from_file("nav1D.seq")
    # sequence3.calculate_stats(stock2)

    sequence4 = Sequence("CypG")
    sequence4.read_from_file("CypG.txt")
    sequence4.calculate_stats(stock4)

    # stock5 = Stock("Wow")
    # stock5.read_from_file("LabelStore_gly_ser_thr.txt")
    # stock5.read_label_prices("prices.txt")
    # stock5.add_usage_dict(USAGE_TABLE)
    #
    # stock6 = Stock("Full")
    # stock6.read_from_file("full_stock.txt")
    # stock6.read_label_prices("prices.txt")
    # stock6.add_usage_dict(USAGE_TABLE)
    #
    # sequence6 = Sequence("Nav2D")
    # sequence6.read_from_file("worst_seq.txt")
    # sequence6.calculate_stats(stock6)
    #
    #
    # sequence5 = Sequence("Nav2D")
    # sequence5.read_from_file("Nav2D.txt")
    # sequence5.calculate_stats(stock2)
    #
    # sequence7 = Sequence("FAM14B")
    # sequence7.read_from_file("FAM14B.seq")
    # sequence7.calculate_stats(stock3)

    solution = Solution("FAM14B", sequence4, stock4)


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
            print(solution)
            #if not solution.found:
                #print(solution)

        if time.time() - t0 > 600 * ten_minute_counter:
            print(str(ten_minute_counter * 10) + " minutes passed")
            print("iteration: " + str(iteration) + "; time elapsed: " + str(round(time.time() - t0, 2)))
            if CHECK_PRICE:
                print("curr best price: " + str(solution.best_price))
            print("")
            ten_minute_counter += 1
            print(solution)

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


def sort_scheme_stats_table(table):
    for i in range(len(table[0])):
        for j in range(len(table[0]) - i -1):
            if int(table[0][j]) > int(table[0][j+1]):
                tmp_0 = table[0][j]
                tmp_1 = table[1][j]
                table[0][j] = table[0][j+1]
                table[1][j] = table[1][j + 1]
                table[0][j+1] = tmp_0
                table[1][j + 1] = tmp_1
    return table


def write_stats_table_to_file(table, filename):
    log_file = open(filename, 'a')
    for i in range(len(table[0])):
        log_file.write(table[0][i] + table[1][i])
    log_file.close()
    return 0


def make_solution_stats(filename):

    stock2 = Stock("Goncharuk")
    stock2.read_from_file("LabelStore_3_ot_seregi.txt")
    stock2.read_label_prices("prices.txt")
    stock2.add_usage_dict(USAGE_TABLE)

    sequence3 = Sequence("Nav1D")
    sequence3.read_from_file("nav1D.seq")
    sequence3.calculate_stats(stock2)

    solution = Solution("Nav1D", sequence3, stock2)
    solution.read_existing_solution("nav1D_best_solution_to_read.txt")
    solution.make_codes_table(filename)


def check_solution(filenames):
    stock = Stock("")
    stock.read_from_file("LabelStore_Maslennikov.txt")
    stock.read_label_prices("prices.txt")

    sequence = Sequence("")
    sequence.read_from_file(filenames[0])
    sequence.calculate_stats(stock)

    solution = Solution("qsec", sequence, stock)

    solution.read_existing_solution(filenames[1])
    solution.print_table()


def clean_sequence(lines):
    joined = "".join(lines)
    sequence = ""
    for char in joined.upper():
        if char in RES_TYPES:
            sequence += char
    return sequence


def get_lines(filename):

    try:
        f = open(filename, 'r', encoding='UTF-8')
        lines = f.readlines()
        f.close()
    except FileNotFoundError:
        return 1
    new_lines = []
    for line in lines:
        curr_line = line.rstrip()
        if curr_line:
            new_lines.append(curr_line)
    return new_lines


def read_config(config_file):

    lines = get_lines(config_file)
    if lines == 1:
        print("ERROR! Config file '" + config_file + "' not found")
        return 1
    task = Task()

    parameters = {}
    for line in lines:
        split_comment = line.split("#")[0]
        split_line = split_comment.split()
        try:
            parameters[split_line[0]] = split_line[1:]
        except IndexError:
            pass

    try:
        job_name = parameters["job_name"]
    except KeyError:
        print("ERROR! Job name is not specified in the config file")
        return 1
    if job_name == []:
        print("ERROR! Job name is not specified in the config file")
        return 1
    task.name = job_name[0]

    try:
        spec_vec = parameters["spectra_types"]
    except KeyError:
        print("ERROR! Spectra types are not specified in the config file")
        return 1
    if spec_vec == []:
        print("ERROR! Spectra types are not specified in the config file")
        return 1
    if task.read_spec_vec(spec_vec) == 1:
        return 1

    try:
        lab_vec = parameters["labeling_types"]
    except KeyError:
        print("ERROR! Labeling types are not specified in the config file")
        return 1
    if lab_vec == []:
        print("ERROR! Labeling types are not specified in the config file")
        return 1
    if task.read_label_vec(lab_vec) == 1:
        return 1

    try:
        stock_file = parameters["stock_file"]
    except KeyError:
        print("ERROR! Stock file is not specified in the config file")
        return 1
    if stock_file == []:
        print("ERROR! Stock file is not specified in the config file")
        return 1
    if task.read_stock(stock_file[0]) == 1:
        return 1

    try:
        optimize_price = parameters["optimize_price"]
    except KeyError:
        print("ERROR! Price optimization is not specified")
        return 1
    if optimize_price[0] == "Yes":
        task.set_price_optimization(True)
        try:
            price_file = parameters["prices_file"]
        except KeyError:
            print("ERROR! Price file is not specified in the config file")
            return 1
        if price_file == []:
            print("ERROR! Price file is not specified in the config file")
            return 1
        if task.read_prices(price_file[0]) == 1:
            return 1
    elif optimize_price[0] == "No":
        task.set_price_optimization(False)
    else:
        print("ERROR! 'optimize_price' field must have value 'Yes' or 'No'")
        return 1

    try:
        sequence_file = parameters["sequence_file"]
    except KeyError:
        print("ERROR! Sequence file is not specified in the config file")
        return 1
    if sequence_file == []:
        print("ERROR! Sequence file is not specified in the config file")
        return 1
    if task.read_sequence(sequence_file[0]) == 1:
        return 1

    return task


def main():
    task = read_config(CONFIG_FILE)
    if task != 1:
        iterations = task.solve()
        if iterations == 1:
            print("Error! Check spectra and labeling types set. All 15N labels are")
        else:
            print("Final number of iterations: ", iterations)
            task.write_results()


if __name__ == '__main__':

    main()
        # while True:
        #     print("Solution is found! Press '1' to exit")
        #     input_key = input()
        #     if input_key == "1":
        #         break