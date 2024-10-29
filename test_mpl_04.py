import matplotlib.pyplot as plt
from matplotlib.widgets import TextBox
import logging
from datetime import datetime

# Configuration du logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class LoggingTextBox:
    def __init__(self, axbox):
        self.text_box = TextBox(axbox, 'Logs', initial='', color='.95', hovercolor='1')
        self.logs = []  # Liste pour stocker tous les logs
        self.display_index = 0  # Index pour la navigation dans les logs
        self.auto_scroll = True  # Drapeau pour le défilement automatique
        # Personnaliser le cadre
        for spine in axbox.spines.values():
            spine.set_edgecolor('black')
            spine.set_linewidth(1.5)

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
        self.text_box.set_val(displayed_logs)
        plt.draw()

    def scroll_logs(self, direction):
        self.auto_scroll = False  # Désactiver le défilement automatique lors de l'utilisation des touches
        if direction == 'up' and self.display_index > 0:
            self.display_index -= 1
        elif direction == 'down' and self.display_index < len(self.logs) - 5:
            self.display_index += 1
        self.display_logs()

# Création de la figure et de l'axe
fig, ax = plt.subplots()
fig.subplots_adjust(bottom=0.3)

# Ajout du TextBox pour les logs
axbox = fig.add_axes([0.1, 0.05, 0.8, 0.2])
log_text_box = LoggingTextBox(axbox)

# Fonction pour ajouter des logs
def log_message(message):
    logger.info(message)
    log_text_box.update_logs(message)

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
        log_text_box.scroll_logs('up')
    elif event.key == 'down':
        log_text_box.scroll_logs('down')

fig.canvas.mpl_connect('key_press_event', on_key_press)

# Exemple de données pour le graphe
x = range(10)
y = [i**2 for i in x]
ax.plot(x, y)

# Affichage de la fenêtre et maintien ouvert
plt.show(block=True)

