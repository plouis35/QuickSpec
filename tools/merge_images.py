import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from astropy.io import fits

# Charger les deux images FITS
img1 = fits.getdata('neon-60s-1.fits')
img2 = fits.getdata('neon-5s-1.fits')

fits_header=fits.getheader('neon-60s-1.fits')

# Vérification des dimensions
assert img1.shape == img2.shape, "Les images doivent avoir la même taille"
#print(fits_header)
height, width = img1.shape

# Fonction de fusion
def fusion_images(x_cut):
    x = int(x_cut)
    fusion = np.zeros_like(img1)
    fusion[:, :x] = img1[:, :x]
    fusion[:, x:] = img2[:, x:]
    return fusion

# Initialisation de l'affichage
fig, (ax_img, ax_prof) = plt.subplots(2, 1, figsize=(10, 6)) #, gridspec_kw={'height_ratios': [3, 1]})
ax_img.axis('off')
fig.suptitle(' ')
plt.subplots_adjust(bottom=0.3)

initial_x = width // 2
fusion = fusion_images(initial_x)

# Affichage image
img_display = ax_img.imshow(fusion, cmap='gray', origin='lower', vmin=np.min(fusion) / 16, vmax=np.max(fusion) / 16, aspect="auto")
ax_img.set_title(f"Fusion à X = {initial_x}")

# Affichage profil spectral
profile_line, = ax_prof.plot(np.sum(fusion, axis=0), color='blue')
#ax_prof.set_yscale('log')
#ax_prof.set_title("Profil spectral (somme par colonne)")
ax_prof.set_xlim(0, width)
ax_prof.set_xlabel("Position X")
ax_prof.set_ylabel("Intensité")

# Slider
ax_slider = plt.axes([0.2, 0.15, 0.6, 0.03])
slider = Slider(ax_slider, 'X', 0, width, valinit=initial_x, valstep=1)

# Bouton SAVE
ax_button = plt.axes([0.45, 0.05, 0.1, 0.04])
button = Button(ax_button, 'Save & exit', color='lightgray', hovercolor='skyblue')

# Mise à jour de l'image et du profil
def update(val):
    x = int(slider.val)
    fusion = fusion_images(x)
    img_display.set_data(fusion)
    ax_img.set_title(f"Fusion à X = {x}")
    profile_line.set_ydata(np.sum(fusion, axis=0))
    fig.canvas.draw_idle()

slider.on_changed(update)

# Sauvegarde de l'image composite
def save(event):
    x = int(slider.val)
    fusion = fusion_images(x)
    fits.writeto('neon_merged.fits', data=fusion, header=fits_header, overwrite=True)
    sys.exit(0)

button.on_clicked(save)

plt.show()
