#Importing the modules and classes  
import tkinter as t  
from tkinter import ttk  
from tkinter.messagebox import showinfo  
  
#Root window creation  
root_window = t.Tk()  
root_window.title('Hierarchical representation')  
root_window.geometry('600x370')  
  
#grid layout  
root_window.rowconfigure(0, weight=1)  
root_window.columnconfigure(0, weight=1)  
  
#Tree view  
tree = ttk.Treeview(root_window)  
tree.heading('#0', text = 'Courses', anchor=t.W)  
  
#Inserting more data  
tree.insert('', t.END, text = 'Engineering', iid=0, open=False)  
tree.insert('', t.END, text = 'Mechanical engineering', iid=1, open=False)  
tree.insert('', t.END, text = 'Electrical engineering', iid=2, open=False)  
tree.insert('', t.END, text = 'Electronic and computer engineering', iid=3, open=False)  
tree.insert('', t.END, text = 'Information Technology', iid=4, open=False)  
tree.insert('', t.END, text = 'Computer science and engineering', iid=5, open=False)  
tree.insert('', t.END, text = 'Artificial Intelligence and Machine Learning', iid=6, open=False)  
tree.insert('', t.END, text = 'Data Science', iid=7, open=False)  
tree.insert('', t.END, text = 'Cyber Security', iid=8, open=False)  
tree.insert('', t.END, text = 'Internet Of Things', iid=9, open=False)  
  
#Assinging the sub-data into Engineering  
tree.move(1, 0, 0)  
tree.move(2, 0, 1)  
tree.move(3, 0, 2)  
tree.move(4, 0, 3)  
tree.move(5, 0, 4)  
  
#Assigning the sub-data to Computer science and engineering  
tree.move(6, 5, 0)  
tree.move(7, 5, 1)  
tree.move(8, 5, 2)  
tree.move(9, 5, 3)  
  
#Positioning  
tree.grid(row = 0, column = 0, sticky = t.NSEW)  
  
#Execution  
root_window.mainloop()  

