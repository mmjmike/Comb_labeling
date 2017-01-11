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

LARGE_FONT = ("Verdana", 12)
NORM_FONT = ("Verdana", 10)
SMALL_FONT = ("Verdana", 8)
STOCK_CELL_SIZE = 10

RES_TYPES = ("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
             "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
LABEL_TYPES = ("X", "N", "C", "D")

label_table = [[False for col in range(20)]
                    for row in range(4)]

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

        self.frames = {}

        for F in (StartPage, MainPage, PricesInput, GraphPage):
            frame = F(container, self)
            self.frames[F] = frame
            frame.grid(row=0, column=0, sticky="nsew")

        self.show_frame(StartPage)

    def default_label(self):
        self.frames[MainPage].label_input.default_label_table()

    def default_prices(self):
        self.frames[PricesInput].table.change(DEFAULT_PRICES_TABLE)

    def show_frame(self, cont):

        frame = self.frames[cont]
        frame.tkraise()


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
        prices_button = ttk.Button(self, text="Input prices",
                            command=lambda: controller.show_frame(PricesInput))
        prices_button.pack()

        label = ttk.Label(self, text='Main page', font=LARGE_FONT)
        label.pack(padx=10, pady=10)
        self.label_input = LabelInput(self)
        self.label_input.pack()
        self.seq_input = SequenceInput(self)
        self.seq_input.pack()

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

        self.text = tk.Text(self, height=6, width=50)
        self.text.pack()

        check_button = ttk.Button(self, text="Use sequence", command=self.check_sequence)
        check_button.pack()
        table_button = ttk.Button(self, text="Show pairs table", command=self.show_table)
        table_button.pack()

    def check_sequence(self):
        global sequence
        sequence = self.text.get("1.0", tk.END)

    def show_table(self):
        pass


class LabelInput(tk.Frame):

    global label_table, edit_labeling_scheme

    def __init__(self, parent):
        tk.Frame.__init__(self, parent)

        self._entry = {}
        self.rows = 4
        self.columns = 20

        self.canvas = tk.Canvas(self, width=600, height=120, bg="white")

        self.edit = tk.BooleanVar()
        self.edit.set(True)

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
        self.canvas.bind("<Button 1>", self._change_table_cell)

        default_button = ttk.Button(self, text="Default stock", command=self.default_label_table)
        default_button.grid(row=6, column=12, columnspan=3)

        read_file_button = ttk.Button(self, text="Read from file",
                                      command=lambda: open_stock_file(parent.parent))
        read_file_button.grid(row=6, column=16, columnspan=3)

        edit_checkbox = ttk.Checkbutton(self, text="Edit",
                                        variable=self.edit, onvalue=True,
                                        offvalue=False)
        edit_checkbox.grid(row=6, column=5, columnspan=4)

    def draw_canvas(self):
        self.canvas.delete("all")

        for i in range(4):
            for j in range(20):
                cell_color = "white"
                if label_table[i][j]:
                    if self.edit.get():
                        cell_color = "green"
                    else:
                        cell_color = "#555555"
                if (i == 1 or i == 3) and j == 12:
                    cell_color = "#999999"
                self.canvas.create_rectangle(j * 30, i * 30, (j + 1) * 30, (i + 1) * 30, fill=cell_color)
        for i in range(3):
            self.canvas.create_line(0, 30 * (i + 1), 600, 30 * (i + 1), width=1, fill="black")
        for i in range(19):
            self.canvas.create_line(30 * (i + 1), 0, 30 * (i + 1), 120, width=1, fill="black")
        self.canvas.after(20, self.draw_canvas)

    def default_label_table(self):
        global label_table

        if self.edit.get():
            for i in range(20):
                for j in range(3):
                    label_table[j][i] = True
            for i in range(20):
                label_table[3][i] = False

    def _change_table_cell(self, event):
        global label_table
        x = event.x // 30
        y = event.y // 30
        if self.edit.get() and x <= 19 and x >= 0 and y <= 3 and y >= 0:
                label_table[y][x] = label_table[y][x] ^ True


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

    app.geometry("800x500")
    app.minsize(650, 430)
    #ani = animation.FuncAnimation(f, animate, interval=1000)
    app.mainloop()
