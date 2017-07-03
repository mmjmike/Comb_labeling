import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import matplotlib.animation as animation
from matplotlib import style
style.use("ggplot")

import tkinter as tk
from tkinter import ttk
import tkinter.constants as Tkconstants
import tkinter.filedialog as tkFileDialog

import file_browser
import csv
import SelLabel

LARGE_FONT = ("Verdana", 12)
NORM_FONT = ("Verdana", 10)
SMALL_FONT = ("Verdana", 8)
STOCK_CELL_SIZE = 10
CELL_HEIGHT = 20
CELL_WIDTH = 50

RES_TYPES = ("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
             "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
LABEL_TYPES = ("X", "N", "C", "D")

label_table = [[False for col in range(20)]
                    for row in range(4)]
label_table[1][12] = False
label_table[3][12] = False
job_name = ""

DEFAULT_PRICES_TABLE = [[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                        [3, 0, 4, 3, 3, 1, 0, 8, 3, 3, 9, 12, 0, 24, 100, 7, 18, 3, 24, 6],
                        [3, 0, 0, 0, 4, 2, 0, 12, 20, 2, 0, 0, 1, 0, 10, 16, 0, 3, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 100, 0, 0, 0, 60, 0, 60, 0, 0, 0, 0, 150, 0]]

prices_table = [[0 for col in range(20)] for row in range(4)]

sequence = ''

START_PAGE_TEXT = """
                   Welcome to SelLabel application!
It is designed to find the optimal solution for combinatorial labeling
       of proteins to facilitate NMR backbone signal assignment.

This version is in process of development. Soon several new options
                will be added for your convenience.

       It is an open source software. Please support us by citing
the following publication if you used this application in your research:

                       ------ publication -------
"""


f = Figure(figsize=(5,5), dpi=100)
a = f.add_subplot(111)


def run(parent):
    run_parameters = {
        "job_name": job_name,
        "sequence": sequence,
        "label_table": label_table,
        "prices_table": prices_table,
        "optimize_price": parent.check_price.get(),
        "HNCA": parent.hnca.get()
    }
    SelLabel.set_parameters(run_parameters)
    SelLabel.main()


def open_stock_file(parent):
    global label_table

    options = {}
    options['defaultextension'] = '.txt'
    options['filetypes'] = [('all files', '.*'), ('text files', '.txt')]
    options['initialdir'] = 'C:\\'
    options['initialfile'] = 'myfile.txt'
    options['parent'] = parent
    options['title'] = 'This is a title'
    label_dict = {}

    filename = tkFileDialog.askopenfilename(**options)
    if filename:
        try:
            with open(filename, 'r') as f:
                reader = csv.reader(f, delimiter='\t')
                d = list(reader)
            # ADD check if file format is correct
            for i in range(len(d) - 1):
                label_dict.update({d[i + 1][0]: ''.join(d[i + 1][1:])})
            for i in range(len(RES_TYPES)):
                for j in range(len(LABEL_TYPES)):
                    if RES_TYPES[i] in label_dict:
                        label_table[j][i] = (label_dict[RES_TYPES[i]][j] == '1')
                    elif j > 0:
                        label_table[j][i] = False
                    else:
                        label_table[j][i] = True
            parent.update_sequence_object()

        except Exception as e:
            popupmsg("Error in reading file: " + str(e))


def open_seq_file(parent):
    global sequence

    options = {}
    options['defaultextension'] = '.seq'
    options['filetypes'] = [('all files', '.*'), ('text files', '.txt')]
    options['initialdir'] = 'C:\\'
    options['initialfile'] = 'myfile.txt'
    options['parent'] = parent
    options['title'] = 'This is a title'

    filename = tkFileDialog.askopenfilename(**options)
    if filename:
        try:
            with open(filename, 'r') as f:
                raw_sequence = f.read()
                f.close()
            if check_sequence(raw_sequence):
                sequence = raw_sequence
                parent.update_sequence_object()
                text_field = parent.frames[MainPage].seq_input.text
                text_field.delete('1.0', tk.END)
                text_field.insert(tk.END, sequence)

        except Exception as e:
            popupmsg("Error in reading file: " + str(e))


def open_prices_file(parent):
    global prices_table

    options = {}
    options['defaultextension'] = '.txt'
    options['filetypes'] = [('all files', '.*'), ('text files', '.txt')]
    options['initialdir'] = 'C:\\'
    options['initialfile'] = 'myfile.txt'
    options['parent'] = parent
    options['title'] = 'This is a title'
    prices_dict = {}

    filename = tkFileDialog.askopenfilename(**options)
    if filename:
        try:
            with open(filename, 'r') as f:
                reader = csv.reader(f, delimiter='\t')
                d = list(reader)
            # ADD check if file format is correct
            for i in range(len(d) - 1):
                residue = d[i + 1][0]
                residue_prices = [
                    d[i + 1][1],
                    d[i + 1][2],
                    d[i + 1][3],
                    d[i + 1][4]
                ]
                prices_dict.update({residue: residue_prices})
            for i in range(len(RES_TYPES)):
                for j in range(len(LABEL_TYPES)):
                    if RES_TYPES[i] in prices_dict:
                        prices_table[j][i] = float(prices_dict[RES_TYPES[i]][j])
                    elif j > 0:
                        prices_table[j][i] = 0
                    else:
                        prices_table[j][i] = 1
            parent.frames[PricesInput].table.change(prices_table)


        except Exception as e:
            popupmsg("Error in reading file: " + str(e))


def check_sequence(sequence):
    return True


def popupmsg(msg):
    popup = tk.Tk()
    popup.wm_title("!")
    label = ttk.Label(popup, text=msg, font=NORM_FONT)
    label.pack(side="top", fill="x", pady=10)
    B1 = ttk.Button(popup, text="OK", command=popup.destroy)
    B1.pack()
    popup.mainloop()


def animate(i):
    pull_data = open("SampleData.txt", "r").read()
    data_list = pull_data.split("\n")
    xList = []
    yList = []
    for each_line in data_list:
        if len(each_line) > 1:
            x, y = each_line.split(',')
            xList.append(int(x))
            yList.append(int(y))

    a.clear()
    a.plot(xList, yList)


class SelLabelapp(tk.Tk):

    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)
        #tk.Tk.iconbitmap(self, default="sellabelicon.ico")
        tk.Tk.wm_title(self, "Sel Label")

        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        menubar = tk.Menu(container)

        filemenu = tk.Menu(menubar, tearoff=0)
        filemenu.add_command(label="Open project",
                             command=lambda: popupmsg("Sorry, this function is not supported yet"))
        filemenu.add_command(label="New project",
                             command=lambda: popupmsg("Sorry, this function is not supported yet"))
        filemenu.add_command(label="Save project",
                             command=lambda: popupmsg("Sorry, this function is not supported yet"))
        filemenu.add_command(label="Save project as...",
                             command=lambda: file_browser.FileBrowser.askopenfile())
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=quit)
        menubar.add_cascade(label="File", menu=filemenu)

        sequence_menu = tk.Menu(menubar, tearoff=0)
        sequence_menu.add_command(label="Read from file...",
                                  command=lambda: open_seq_file(self))
        sequence_menu.add_command(label="Save as...",
                                  command=lambda: popupmsg("Sorry, this function is not supported yet"))
        menubar.add_cascade(label="Sequence", menu=sequence_menu)

        stock_menu = tk.Menu(menubar, tearoff=0)
        stock_menu.add_command(label="Read from file...",
                               command=lambda: open_stock_file(self))
        stock_menu.add_command(label="Read prices from file...",
                               command=lambda: open_prices_file(self))
        stock_menu.add_separator()
        stock_menu.add_command(label="Default stock",
                               command=self.default_label)
        stock_menu.add_command(label="Default prices",
                               command=self.default_prices)
        stock_menu.add_separator()
        stock_menu.add_command(label="Save as...",
                               command=lambda: popupmsg("Sorry, this function is not supported yet"))
        stock_menu.add_command(label="Save prices as...",
                               command=lambda: popupmsg("Sorry, this function is not supported yet"))
        menubar.add_cascade(label="Stock", menu=stock_menu)

        tk.Tk.config(self, menu=menubar)

        self.panels = (StartPage, MainPage, PricesInput, GraphPage, AllPairsTable, StockPairsTable
                       # , LabelInput
                       )

        self.frames = {}

        self.first_run = True

        self.update_sequence_object()

        for F in self.panels:
            self.frame = F(container, self)
            self.frames[F] = self.frame
            self.frame.grid(row=2, column=0, sticky="nsew")

        self.show_frame(MainPage)

    def default_label(self):
        self.frames[MainPage].label_input.default_label_table()

    def default_prices(self):
        self.frames[PricesInput].table.change(DEFAULT_PRICES_TABLE)

    def show_frame(self, cont):

        if self.frame != self.frames[cont]:
            self.frame = self.frames[cont]
            self.frame.tkraise()

    def update_sequence_object(self):
        global sequence

        self.sequence_obj = SelLabel.Sequence("Sequence", sequence)
        stock_obj = SelLabel.Stock("stock")
        stock_obj.read_from_table(label_table)
        self.sequence_obj.calculate_stats(stock_obj)
        if not self.first_run:
            self.frames[AllPairsTable].draw_canvas()
            self.frames[StockPairsTable].draw_canvas()
        self.first_run = False


def qf(param):
    print(param)


class StartPage(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = ttk.Label(self, text=START_PAGE_TEXT, font=LARGE_FONT)
        label.pack(padx=10, pady=10)

        agree_button = ttk.Button(self, text="I Agree",
                           command=lambda : controller.show_frame(MainPage))
        agree_button.pack()

        exit_button = ttk.Button(self, text="Exit", command=quit)
        exit_button.pack()


class MainPage(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.parent = parent
        self.controller = controller

        button_panel = ButtonPanel(self)
        button_panel.pack(fill=tk.X)
        ttk.Separator(self).pack(fill=tk.X)

        ttk.Label(self, text='Sequence', font=LARGE_FONT).pack()
        # self.label_input = LabelInput(self)
        # self.label_input.grid(row=3, column=0, columnspan=4, sticky="ew")
        self.seq_input = SequenceInput(self)
        self.seq_input.pack()

        self.edit = tk.BooleanVar()
        self.edit.set(True)

        self.edit_select = tk.IntVar()

        tk.Radiobutton(self, text="Edit stock", variable=self.edit_select, value=1).pack()
        tk.Radiobutton(self, text="Edit prices", variable=self.edit_select, value=2).pack()
        tk.Radiobutton(self, text="Lock edit", variable=self.edit_select, value=3).pack()
        self.edit_select.set(1)

        self.label_canvas = LabelCanvas(self)
        self.label_canvas.pack()

        default_button = ttk.Button(self, text="Default stock", command=self.label_canvas.default_label_table)
        default_button.pack()

        read_file_button = ttk.Button(self, text="Read from file",
                                      command=lambda: open_stock_file(controller))
        read_file_button.pack()

        edit_checkbox = ttk.Checkbutton(self, text="Edit",
                                        variable=self.edit, onvalue=True,
                                        offvalue=False)
        edit_checkbox.pack()


        self.check_price = tk.BooleanVar()
        self.check_price.set(False)
        price_checkbox = ttk.Checkbutton(self, text="Optimize price",
                                         variable=self.check_price, onvalue=True, offvalue=False)
        price_checkbox.pack()

        self.hnca = tk.BooleanVar()
        self.hnca.set(False)
        hnca_checkbox = ttk.Checkbutton(self, text="Use HNCA",
                                        variable=self.hnca, onvalue=True, offvalue=False)
        hnca_checkbox.pack()

        run_button = ttk.Button(self, text="Run", command=lambda: run(self))
        run_button.pack()


class SequenceInput(tk.Frame):

    def __init__(self, parent):
        tk.Frame.__init__(self, parent)
        self.parent = parent

        self.text = tk.Text(self, height=6, width=50)
        self.text.pack()

        check_button = ttk.Button(self, text="Use sequence", command=self.check_sequence)
        check_button.pack()

    def check_sequence(self):
        global sequence
        sequence = self.text.get("1.0", tk.END)[:-1].upper()
        self.parent.controller.update_sequence_object()


    def show_table(self):
        global sequence
        test_sequence = self.text.get("1.0", tk.END).upper()
        test_sequence = test_sequence[:-1]

        for i in range(len(test_sequence)):
            if test_sequence[i] not in RES_TYPES and test_sequence[i] != "-":
                message = "Error! Wrong character in sequence (position {0})".format(i)
                popupmsg(message)
                return
        pair_table()


def pair_table():
    p_table = PairTable()
    p_table.wm_title("Table of pairs")
    p_table.mainloop()


class PairTable(tk.Tk):

    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)
        tk.Tk.wm_title(self, "Pair Table")

        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.frames = {}

        for F in (AllPairsTable, StockPairsTable):
            frame = F(container, self)
            self.frames[F] = frame
            frame.grid(row=0, column=0, sticky="nsew")

        self.show_frame(AllPairsTable)

    def show_frame(self, cont):
        frame = self.frames[cont]
        frame.tkraise()


class AllPairsTable(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        button_panel = ButtonPanel(self)
        button_panel.pack(fill=tk.X)
        ttk.Separator(self).pack(fill=tk.X)

        label = ttk.Label(self, text='All pairs in sequence sorted', font=LARGE_FONT)
        label.pack()

        self.canvas = tk.Canvas(self, width=CELL_WIDTH, height=CELL_HEIGHT,
                                bg="black")
        self.canvas.pack()

        self.stats_label = ttk.Label(self, text='', font=LARGE_FONT)
        self.stats_label.pack()

        self.draw_canvas()

    def draw_canvas(self):
        if self.controller.sequence_obj.sequence != "":
            row_res = self.controller.sequence_obj.residues_first
            col_res = self.controller.sequence_obj.residues_second
            x = len(col_res)
            y = len(row_res)
            cells = (x, y)
            self.canvas.config(width=CELL_WIDTH * (cells[0] + 1), height=CELL_HEIGHT * (cells[1] + 1))
            res_pairs = self.controller.sequence_obj.all_residue_pairs
            unique_pairs_count = 0
            pairs_count = 0

            self.canvas.delete("all")
            for i in range(cells[1]):
                self.canvas.create_line(0, CELL_HEIGHT * (i + 1), CELL_WIDTH * (cells[0] + 1), CELL_HEIGHT * (i + 1),
                              fill="white")
                self.canvas.create_text(CELL_WIDTH / 2, CELL_HEIGHT * (i + 1.5), text=row_res[i], fill="white")
            for i in range(cells[0]):
                self.canvas.create_line(CELL_WIDTH * (i + 1), 0, CELL_WIDTH * (i + 1), CELL_HEIGHT * (cells[1] + 1),
                              fill="white")
                self.canvas.create_text(CELL_WIDTH * (i + 1.5), CELL_HEIGHT / 2, text=col_res[i], fill="white")
            for i in range(cells[1]):
                for j in range(cells[0]):
                    if res_pairs[i][j] > 0:
                        if res_pairs[i][j] == 1:
                            color = "green"
                            unique_pairs_count += 1
                        else:
                            pairs_count += 1
                            color = "red"
                        self.canvas.create_rectangle(CELL_WIDTH * (j + 1) + 1,
                                                     CELL_HEIGHT * (i + 1) + 1,
                                                     CELL_WIDTH * (j + 2) - 1,
                                                     CELL_HEIGHT * (i + 2) - 1,
                                                     fill=color)
                        self.canvas.create_text(CELL_WIDTH * (j + 1.5),
                                                CELL_HEIGHT * (i + 1.5),
                                                text=res_pairs[i][j],
                                                fill="white")
            all_pairs_count = unique_pairs_count + pairs_count
            stats_text = "Unique: {}, Non-unique: {}, Total: {}".format(str(unique_pairs_count),
                                                                         str(pairs_count), str(all_pairs_count))
            self.stats_label.config(text=stats_text)
        else:
            self.stats_label.config(text="")
            self.canvas.config(width=CELL_WIDTH, height=CELL_HEIGHT)


class StockPairsTable(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller

        button_panel = ButtonPanel(self)
        button_panel.pack(fill=tk.X)
        ttk.Separator(self).pack(fill=tk.X)

        label = ttk.Label(self, text='Prices Input', font=LARGE_FONT)
        label.pack()

        self.canvas = tk.Canvas(self, width=CELL_WIDTH, height=CELL_HEIGHT,
                                bg="black")
        self.canvas.pack()

        self.stats_label = ttk.Label(self, text='', font=LARGE_FONT)
        self.stats_label.pack()

        self.other_label = ttk.Label(self, text='', font=LARGE_FONT)
        self.other_label.pack()

        self.draw_canvas()

    def draw_canvas(self):
        if self.controller.sequence_obj.sequence != "":
            row_res = self.controller.sequence_obj.residues_carbon
            col_res = self.controller.sequence_obj.residues_nitro
            x = len(col_res)
            y = len(row_res)
            cells = (x, y)
            self.canvas.config(width=CELL_WIDTH * (cells[0] + 1), height=CELL_HEIGHT * (cells[1] + 1))
            res_pairs = self.controller.sequence_obj.residue_pairs
            unique_pairs_count = 0
            pairs_count = 0
            other_nitro_list = self.controller.sequence_obj.residues_not_nitro
            other_carbon_list = self.controller.sequence_obj.residues_not_carbon

            self.canvas.delete("all")
            for i in range(cells[1]):
                self.canvas.create_line(0, CELL_HEIGHT * (i + 1), CELL_WIDTH * (cells[0] + 1), CELL_HEIGHT * (i + 1),
                              fill="white")
                self.canvas.create_text(CELL_WIDTH / 2, CELL_HEIGHT * (i + 1.5), text=row_res[i], fill="white")
            for i in range(cells[0]):
                self.canvas.create_line(CELL_WIDTH * (i + 1), 0, CELL_WIDTH * (i + 1), CELL_HEIGHT * (cells[1] + 1),
                              fill="white")
                self.canvas.create_text(CELL_WIDTH * (i + 1.5), CELL_HEIGHT / 2, text=col_res[i], fill="white")
            for i in range(cells[1]):
                for j in range(cells[0]):
                    if res_pairs[i][j] > 0:
                        if res_pairs[i][j] == 1:
                            color = "green"
                            unique_pairs_count += 1
                        else:
                            color = "red"
                            pairs_count += 1
                        self.canvas.create_rectangle(CELL_WIDTH * (j + 1) + 1,
                                                     CELL_HEIGHT * (i + 1) + 1,
                                                     CELL_WIDTH * (j + 2) - 1,
                                                     CELL_HEIGHT * (i + 2) - 1,
                                                     fill=color)
                        self.canvas.create_text(CELL_WIDTH * (j + 1.5),
                                                CELL_HEIGHT * (i + 1.5),
                                                text=res_pairs[i][j],
                                                fill="white")
            all_pairs_count = unique_pairs_count + pairs_count
            stats_text = "Unique: {}, Non-unique: {}, Total: {}".format(str(unique_pairs_count),
                                                                         str(pairs_count), str(all_pairs_count))
            other_text = "Other carbon: {}; Other nitro: {}".format(",".join(other_carbon_list),
                                                                    ",".join(other_nitro_list))
            self.stats_label.config(text=stats_text)
            self.other_label.config(text=other_text)
        else:
            self.stats_label.config(text="")
            self.other_label.config(text="")
            self.canvas.config(width=CELL_WIDTH, height=CELL_HEIGHT)


class LabelInput(tk.Frame):

    global label_table, edit_labeling_scheme

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        self.parent = parent

        button_panel = ButtonPanel(self)
        button_panel.pack(fill=tk.X)
        ttk.Separator(self).pack(fill=tk.X)

        self.edit = tk.BooleanVar()
        self.edit.set(True)

        self.label_canvas = LabelCanvas(self)
        self.label_canvas.pack()

        default_button = ttk.Button(self, text="Default stock", command=self.label_canvas.default_label_table)
        default_button.pack()

        read_file_button = ttk.Button(self, text="Read from file",
                                      command=lambda: open_stock_file(controller))
        read_file_button.pack()

        edit_checkbox = ttk.Checkbutton(self, text="Edit",
                                        variable=self.edit, onvalue=True,
                                        offvalue=False)
        edit_checkbox.pack()


class LabelCanvas(tk.Frame):

    def __init__(self, parent):
        tk.Frame.__init__(self, parent)
        self.parent = parent
        self.start = True
        self.rows = 4
        self.columns = 20
        self.stock_cell_width = 45
        self.stock_cell_height = 30
        self.canvas = tk.Canvas(self, width=self.stock_cell_width*self.columns,
                                height=self.stock_cell_height*self.rows, bg="white")
        self.canvas.bind("<Button 1>", self._change_table_cell)
        self.default_label_table()
        self.draw_canvas()



        for row in range(self.rows):
            ttk.Label(self, text=LABEL_TYPES[row], font=NORM_FONT).grid(row=(row + 1), column=0)
        for column in range(self.columns):
            ttk.Label(self, text=RES_TYPES[column], font=NORM_FONT).grid(row=0, column=(column + 1))
        self.canvas.grid(row=1, column=1, columnspan=20, rowspan=4)
        for column in range(self.columns):
            self.grid_columnconfigure(column + 1, weight=1)
        # designate a final, empty row to fill up any extra space
        self.grid_rowconfigure(self.rows + 1, weight=1)

    def draw_canvas(self):

        self.canvas.delete("all")

        for i in range(4):
            for j in range(20):
                cell_color = "white"
                if label_table[i][j]:
                    if self.parent.edit_select.get() != 3:
                        cell_color = "green"
                    else:
                        cell_color = "#555555"
                if (i == 1 or i == 3) and j == 12:
                    cell_color = "#999999"
                self.canvas.create_rectangle(j * self.stock_cell_width,
                                             i * self.stock_cell_height,
                                             (j + 1) * self.stock_cell_width,
                                             (i + 1) * self.stock_cell_height,
                                             fill=cell_color)
        for i in range(3):
            self.canvas.create_line(0, self.stock_cell_height * (i + 1),
                                    self.stock_cell_width * self.columns, self.stock_cell_height * (i + 1), width=1,
                                    fill="black")
        for i in range(19):
            self.canvas.create_line(self.stock_cell_width * (i + 1), 0, self.stock_cell_width * (i + 1),
                                    self.stock_cell_height * self.rows, width=1, fill="black")
        self.canvas.after(20, self.draw_canvas)
        self.canvas.create_text(20, 20, text="1578.83", fill="black")

    def default_label_table(self):
        global label_table

        if self.parent.edit.get():
            for i in range(20):
                for j in range(3):
                    label_table[j][i] = True
            for i in range(20):
                label_table[3][i] = False
        if not self.start:
            self.parent.controller.update_sequence_object()
            self.start = False

    def _change_table_cell(self, event):
        global label_table
        x = event.x // self.stock_cell_width
        y = event.y // self.stock_cell_height
        if self.parent.edit_select.get() == 1 and x <= 19 and x >= 0 and y <= 3 and y >= 0:
            label_table[y][x] = label_table[y][x] ^ True
        if self.parent.edit_select.get() == 2 and x <= 19 and x >= 0 and y <= 3 and y >= 0:
            pass
        label_table[1][12] = False
        label_table[3][12] = False
        self.parent.controller.update_sequence_object()


class PricesInput(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        back_button = ttk.Button(self, text="Back",
                           command=lambda: controller.show_frame(MainPage))
        back_button.pack()

        label = ttk.Label(self, text='Prices Input', font=LARGE_FONT)
        label.pack(padx=10, pady=10)

        self.table = TableInput(self)
        self.table.change(prices_table)
        self.submit = ttk.Button(self, text="Submit", command=self.on_submit)
        self.default = ttk.Button(self, text="Default", command=self.on_default)
        self.zero = ttk.Button(self, text="Clear", command=self.on_zero)
        self.table.pack(side="top", fill="both", expand=True)
        self.submit.pack(side="bottom")
        self.default.pack(side="bottom")
        self.zero.pack(side="bottom")

    def on_submit(self):
        global prices_table
        prices_table = self.table.get()

    def on_default(self):
        global prices_table
        prices_table = DEFAULT_PRICES_TABLE
        self.table.change(prices_table)

    def on_zero(self):
        global prices_table
        prices_table = [[0 for col in range(20)]
                    for row in range(4)]
        self.table.change(prices_table)


class TableInput(tk.Frame):

    def __init__(self, parent):
        tk.Frame.__init__(self, parent)

        self._entry = {}
        self.rows = 4
        self.columns = 20

        # register a command to use for validation
        vcmd = (self.register(self._validate), "%P")

        # create table header
        for column in range(self.columns):
            ttk.Label(self, text=RES_TYPES[column], font=NORM_FONT).grid(row=0, column=(column+1))

        # create the table of widgets
        for row in range(self.rows):
            ttk.Label(self, text=LABEL_TYPES[row], font=NORM_FONT).grid(row=(row+1), column=0)
            for column in range(self.columns):
                index = (row, column)
                e = ttk.Entry(self, validate="key", validatecommand=vcmd, justify=tk.RIGHT)
                e.grid(row=(row+1), column=(column+1))
                self._entry[index] = e
        # adjust column weights so they all expand equally
        for column in range(self.columns):
            self.grid_columnconfigure(column + 1, weight=1)
        # designate a final, empty row to fill up any extra space
        self.grid_rowconfigure(self.rows + 1, weight=1)

    def get(self):
        '''Return a list of lists, containing the data in the table'''
        result = []
        for row in range(self.rows):
            current_row = []
            for column in range(self.columns):
                index = (row, column)
                current_row.append(float(self._entry[index].get()))
            result.append(current_row)
        return result

    def change(self, table):
        for row in range(self.rows):
            for column in range(self.columns):
                index = (row, column)
                self._entry[index].delete(0, tk.END)
                self._entry[index].insert(0, str(float(table[row][column])))

    def _validate(self, P):
        '''Perform input validation.

        Allow only an empty value, or a value that can be converted to a float
        '''
        if P.strip() == "":
            return True

        try:
            f = float(P)
        except ValueError:
            self.bell()
            return False
        return True


class ButtonPanel(tk.Frame):

    def __init__(self, parent):
        tk.Frame.__init__(self, parent)

        self.styles = {}
        for panel in parent.controller.panels:
            self.styles[panel] = "TButton"
        self.styles[parent.__class__] = "current.TButton"

        ttk.Style().configure("current.TButton", background="black")

        sequence_button = ttk.Button(self, text="Main page", style=self.styles[MainPage],
                                     command=lambda: parent.controller.show_frame(MainPage))
        sequence_button.grid(row=0, column=0, sticky='w')
        # stock_button = ttk.Button(self, text="Stock", style=self.styles[LabelInput],
        #                           command=lambda: parent.controller.show_frame(LabelInput))
        # stock_button.grid(row=0, column=1, sticky='w')
        all_paires_button = ttk.Button(self, text="All pairs table", style=self.styles[AllPairsTable],
                                       command=lambda: parent.controller.show_frame(AllPairsTable))
        all_paires_button.grid(row=0, column=2, sticky='w')
        paires_button = ttk.Button(self, text="Stock pairs table", style=self.styles[StockPairsTable],
                                   command=lambda: parent.controller.show_frame(StockPairsTable))
        paires_button.grid(row=0, column=3, sticky='w')


class GraphPage(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = ttk.Label(self, text='Graph Page', font=LARGE_FONT)
        label.pack(padx=10, pady=10)

        button = ttk.Button(self, text="Go to Start Page",
                           command=lambda: controller.show_frame(StartPage))
        button.pack()

        canvas = FigureCanvasTkAgg(f, self)
        canvas.show()

        toolbar = NavigationToolbar2TkAgg(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)


if __name__ == '__main__':
    app = SelLabelapp()

    app.geometry("1000x500")
    app.minsize(650, 430)
    #ani = animation.FuncAnimation(f, animate, interval=1000)
    app.mainloop()
