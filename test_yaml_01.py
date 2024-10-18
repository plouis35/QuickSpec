import tkinter as tk
from tkinter import ttk, filedialog, simpledialog
import yaml

class YAMLTreeView(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("YAML TreeView")
        self.geometry("800x600")

        self.tree = ttk.Treeview(self)
        self.tree.pack(expand=1, fill='both')

        self.tree.bind('<Double-1>', self.on_double_click)

        self.menu = tk.Menu(self)
        self.config(menu=self.menu)
        file_menu = tk.Menu(self.menu, tearoff=0)
        self.menu.add_cascade(label="Fichier", menu=file_menu)
        file_menu.add_command(label="Ouvrir", command=self.open_file)
        file_menu.add_command(label="Quitter", command=self.quit)

    def open_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("YAML files", "*.yaml")])
        if file_path:
            with open(file_path, 'r') as file:
                content = yaml.safe_load(file)
                self.populate_tree(content)

    def populate_tree(self, content, parent=''):
        for key, value in content.items():
            node = self.tree.insert(parent, 'end', text=key)
            if isinstance(value, dict):
                self.populate_tree(value, node)
            else:
                self.tree.insert(node, 'end', text=value)

    def on_double_click(self, event):
        item = self.tree.selection()[0]
        parent = self.tree.parent(item)
        if parent:  # Only allow editing of the value nodes
            current_value = self.tree.item(item, 'text')
            new_value = simpledialog.askstring("Modifier la valeur", "Nouvelle valeur:", initialvalue=current_value)
            if new_value is not None:
                self.tree.item(item, text=new_value)

if __name__ == "__main__":
    app = YAMLTreeView()
    app.mainloop()
