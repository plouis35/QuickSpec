#Importing the modules and the classes  
from tkinter import ttk   
import tkinter as t  
  
#Creating a window   
root_window = t.Tk()  
root_window.resizable(width = 2, height = 2)  
  
#Treeview widget  
tree = ttk.Treeview(root_window, selectmode = 'browse')  
tree.pack(side = 'left')  
  
#scrollbar  
scrollbar = ttk.Scrollbar(root_window, orient = "vertical", command = tree.yview)  
scrollbar.pack(side = 'left', fill = 'x')  
tree.configure(xscrollcommand = scrollbar.set)  
  
#Number of columns   
tree["columns"] = ("1", "2")  
   
# Headings  
tree['show'] = 'headings'  
   
#Columns configuration  
tree.column("1", width = 90, anchor ='c')  
tree.column("2", width = 90, anchor ='se')  
  
#Headings  
tree.heading("1", text = "Name")  
tree.heading("2", text = "Age")  
  
#Data  
tree.insert("", 'end', text ="L1",  
             values =("Harry Styles", "28"))  
tree.insert("", 'end', text ="L2",  
             values =("Zayn Malik", "29"))  
tree.insert("", 'end', text ="L3",  
             values =("Liam Payne", "29"))  
tree.insert("", 'end', text ="L4",  
             values =("Louis Tomlinson", "30"))  
tree.insert("", 'end', text ="L5",  
             values =("Niall Horan", "29"))  
tree.insert("", 'end', text ="L6",  
             values =("Adam Levine", "43"))  
tree.insert("", 'end', text ="L7",  
             values =("Selena Gomez", "30"))  
tree.insert("", 'end', text ="L8",  
             values =("Shawn Mendes", "24"))  
tree.insert("", 'end', text ="L10",  
             values =("Taylor Swift", "32"))  
tree.insert("", 'end', text ="L11",  
             values =("Justin Beiber", "28"))  
tree.insert("", 'end', text ="L12",  
             values =("Demi Lovato", "30"))  
tree.insert("", 'end', text ="L13",  
             values =("Conan Gray", "23"))  
   
# Calling main loop  
root_window.mainloop()  


