import tkinter as tk
from tkinter import scrolledtext
import matplotlib.pyplot as plt
import logging
from datetime import datetime

# Configuration du logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class LoggingWindow:
    def __init__(self, root):
        self.text_box = scrolledtext.ScrolledText(root, wrap=tk.WORD, state='disabled', width=80, height=20)
        self.text_box.pack(padx=10, pady=10)
        self.logs = []  # Liste pour stocker tous les logs
        self.display_index = 0  # Index pour la navigation dans les logs
        self.auto_scroll = True  # Drapeau pour le défilement automatique

    def update_logs(self, message):
        current_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        log_entry = f'{current_time} - {message}'
        self.logs.append(log_entry)
        if self.auto_scroll:
            self.display_index = max(0, len(self.logs) - 5)
        self.display_logs()

    def display_logs(self):
        # Limiter l'affichage aux 5 dernières lignes visibles
        displayed_logs = '\n'.join(self.logs[self.display_index:self.display_index + 5])
        self.text_box.config(state='normal')
        self.text_box.delete(1.0, tk.END)
        self.text_box.insert(tk.END, displayed_logs)
        self.text_box.config(state='disabled')
        self.text_box.yview(tk.END)  # Scroll to the end

    def scroll_logs(self, direction):
        self.auto_scroll = False  # Désactiver le défilement automatique lors de l'utilisation des touches
        if direction == 'up' and self.display_index > 0:
            self.display_index -= 1
        elif direction == 'down' and self.display_index < len(self.logs) - 5:
            self.display_index += 1
        self.display_logs()

# Création de la fenêtre principale pour les logs
root = tk.Tk()
root.title("Logging Window")
log_window = LoggingWindow(root)

# Création de la figure et de l'axe pour le graphique
fig, ax = plt.subplots()
fig.subplots_adjust(bottom=0.3)

# Fonction pour ajouter des logs
def log_message(message):
    logger.info(message)
    log_window.update_logs(message)

# Exemple d'utilisation
log_message("Démarrage du programme...")
log_message("Chargement des données...")
log_message("Traitement en cours...")

# Ajout d'une ligne verticale interactive
vertical_line = ax.axvline(x=0, color='gray', linestyle='--')
x_text = ax.text(0, 1, '', transform=ax.get_xaxis_transform(), ha='left', va='bottom')
y_text = ax.text(0, 0, '', transform=ax.get_yaxis_transform(), ha='left', va='bottom')

def on_mouse_move(event):
    if event.inaxes == ax:
        vertical_line.set_xdata([event.xdata])
        x_text.set_position((event.xdata, 1))
        y_text.set_position((event.xdata, event.ydata))
        x_text.set_text(f'X: {event.xdata:.2f}')
        y_text.set_text(f'Y: {event.ydata:.2f}')
        ax.figure.canvas.draw()
        log_message(f'X: {event.xdata:.2f}, Y: {event.ydata:.2f}')

fig.canvas.mpl_connect('motion_notify_event', on_mouse_move)

# Gestion des touches pour le défilement des logs
def on_key_press(event):
    if event.key == 'up':
        log_window.scroll_logs('up')
    elif event.key == 'down':
        log_window.scroll_logs('down')

fig.canvas.mpl_connect('key_press_event', on_key_press)

# Exemple de données pour le graphe
x = range(10)
y = [i**2 for i in x]
ax.plot(x, y)

# Affichage de la fenêtre de logs et du graphique
plt.show(block=False)
root.mainloop()

