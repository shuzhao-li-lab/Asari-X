from default_parameters import PARAMETERS
import main
import tkinter as tk
from tkinter import scrolledtext, messagebox, filedialog
import threading
from contextlib import redirect_stdout, redirect_stderr

# todo - refactor to use only a single tk root.
def select_directory(params):
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(title="Select Input Directory")
    params['input']['value'] = directory
    root.destroy()
    return params

def parameter_selection(params):
    result = {}
    def run_callback():
        nonlocal result
        output = {}
        for key, widget in widgets.items():
            if types[key] is bool:
                output[key] = widget.get()
            else:
                val = widget.get()
                if types[key] is int:
                    try:
                        output[key] = int(val)
                    except ValueError:
                        messagebox.showerror("Error", f"Invalid integer for {key}")
                        return
                elif types[key] is float:
                    try:
                        output[key] = float(val)
                    except ValueError:
                        messagebox.showerror("Error", f"Invalid float for {key}")
                        return
                else:
                    output[key] = val
        result = output
        root.destroy()

    root = tk.Tk()
    root.title("Edit Parameters")

    canvas = tk.Canvas(root)
    vsb = tk.Scrollbar(root, orient="vertical", command=canvas.yview)
    canvas.configure(yscrollcommand=vsb.set)
    vsb.pack(side="right", fill="y")
    canvas.pack(side="left", fill="both", expand=True)

    frame = tk.Frame(canvas)
    canvas.create_window((0, 0), window=frame, anchor="nw")
    frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))

    widgets = {}
    types = {}
    row = 0
    for key, value in params.items():
        tk.Label(frame, text=key).grid(row=row, column=0, padx=5, pady=2, sticky="w")
        try:
            to_display = value['default'] if 'value' not in value else value['value']
            if value["types"][0] in [bool]:
                var = tk.BooleanVar(value=to_display)
                chk = tk.Checkbutton(frame, variable=var)
                chk.grid(row=row, column=1, padx=5, pady=2, sticky="w")
                widgets[key] = var
                types[key] = bool
            elif value["types"][0] in [int, float]:
                entry = tk.Entry(frame)
                entry.insert(0, str(to_display))
                entry.grid(row=row, column=1, padx=5, pady=2, sticky="w")
                widgets[key] = entry
                types[key] = type(to_display)
            elif value["types"][0] in [str, type(None)]:
                entry = tk.Entry(frame)
                if to_display is None:
                    entry.insert(0, "NONE")
                else:
                    entry.insert(0, to_display)
                entry.grid(row=row, column=1, padx=5, pady=2, sticky="w")
                widgets[key] = entry
                types[key] = str
            row += 1
        except:
            pass

    run_button = tk.Button(root, text="Continue", command=run_callback)
    run_button.pack(side="bottom", pady=10)
    root.mainloop()
    for k, v in result.items():
        if v == "NONE":
            result[k] = None
    return result

def subcommand_selection(params):
    result = {}
    def set_value(value):
        nonlocal result
        params['run'] = value
        result = params
        root.destroy()

    root = tk.Tk()
    root.title("Select an Option")
    for option in main.main({}, dry_run=True):
        btn = tk.Button(root, text=option, command=lambda opt=option: set_value(opt))
        btn.pack(pady=2)
    root.mainloop()
    return result
    
def gui_main(params):
    params = select_directory(params)
    params = parameter_selection(params)
    params = subcommand_selection(params)
    
if __name__ == '__main__':
    gui_main(PARAMETERS)