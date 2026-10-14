import logging

# Constante physique pour la bande Tellurique O2 (Bande A)
LAMBDA_O2_NM = 760.0  # 7600 Å (bande d'absorption atmosphérique de l'oxygène)

def compute_sa_distance_from_o2(num_pixels: float,
                                pixel_size_um: float = 4.8,
                                binning: int = 1,
                                lines_per_mm: float = 200.0) -> float:
    """
    Calcule la distance SA -> capteur (mm) depuis la bande O2 à 7600 Å.

    Args:
        num_pixels (float): Nombre de pixels entre l'Ordre 0 et le creux O2
        pixel_size_um (float): Taille native du pixel en µm (ex: 4.8)
        binning (int): Facteur de binning (1 pour 1x1, 2 pour 2x2, etc.)
        lines_per_mm (float): Traits/mm du réseau (200.0 pour SA200)

    Returns:
        float: Distance SA -> Capteur en mm
    """
    # Taille effective du pixel binné
    effective_pixel_um = pixel_size_um * binning

    # Déplacement physique sur le capteur (mm)
    x_sensor_mm = (num_pixels * effective_pixel_um) / 1000.0

    # Pas du réseau d en mm (ex: 1/200 = 0.005 mm)
    d_mm = 1.0 / lines_per_mm

    # Convertit 760 nm en mm (760e-6 mm)
    lambda_o2_mm = LAMBDA_O2_NM * 1e-6

    # Distance D = (d * x) / λ
    distance_mm = (d_mm * x_sensor_mm) / lambda_o2_mm

    dispersion_AA_per_px = (LAMBDA_O2_NM * 10.0) / num_pixels

    print(f"--- RÉSULTATS DU CALCUL (Bande O2 @ 7600 Å) ---")
    print(f"Pixels Ordre 0 -> O2 (7600 Å)    : {num_pixels:.2f} px")
    print(f"Taille de pixel effective          : {effective_pixel_um:.2f} µm (Binning {binning}x{binning})")
    print(f"Dispersion mesurée                : {dispersion_AA_per_px:.2f} Å/px")
    print(f"-----------------------------------------------")
    print(f"Distance SA200 -> Capteur         : {distance_mm:.2f} mm")

    return distance_mm


# Constantes physiques
LAMBDA_HALPHA_NM = 656.28  # Longueur d'onde H-alpha (6562.8 Å)

def compute_sa_distance_mm(num_pixels: float, 
                           pixel_size_um: float = 4.8, 
                           binning: int = 1, 
                           lines_per_mm: float = 200.0) -> float:
    """
    Calcule la distance réseau SA -> capteur en mm à partir de la raie H-alpha.

    Args:
        num_pixels (float): Nombre de pixels entre l'Ordre 0 et H-alpha
        pixel_size_um (float): Taille native du pixel en µm (ex: 4.8)
        binning (int): Facteur de binning (1 pour 1x1, 2 pour 2x2, etc.)
        lines_per_mm (float): Nombre de traits/mm du réseau (200.0 pour SA200)

    Returns:
        float: Distance SA -> capteur en mm
    """
    # Taille effective du pixel sur l'image traitée
    effective_pixel_um = pixel_size_um * binning

    # Distance physique sur le capteur entre l'Ordre 0 et H-alpha (en mm)
    x_sensor_mm = (num_pixels * effective_pixel_um) / 1000.0

    # Pas du réseau d en mm (ex: 1/200 = 0.005 mm)
    d_mm = 1.0 / lines_per_mm

    # Loi du réseau en approximation linéaire (faibles angles) :
    # sin(theta) ≈ tan(theta) = x / D  =>  λ = d * (x / D)  =>  D = (d * x) / λ
    # Avec λ_Halpha en mm : 656.28 nm = 656.28e-6 mm
    lambda_halpha_mm = LAMBDA_HALPHA_NM * 1e-6

    distance_mm = (d_mm * x_sensor_mm) / lambda_halpha_mm

    # Informations complémentaires
    dispersion_AA_per_px = (LAMBDA_HALPHA_NM * 10.0) / num_pixels

    print(f"--- RÉSULTATS DU CALCUL ---")
    print(f"Distance mesurée Ordre 0 -> H-alpha : {num_pixels:.2f} px")
    print(f"Taille de pixel effective          : {effective_pixel_um:.2f} µm (Binning {binning}x{binning})")
    print(f"Dispersion mesurée sur le spectre  : {dispersion_AA_per_px:.2f} Å/px")
    print(f"-----------------------------------")
    print(f"Distance SA200 -> Capteur         : {distance_mm:.2f} mm")

    return distance_mm

# --- EXEMPLE D'UTILISATION ---
if __name__ == "__main__":
    # Remplace cette valeur par ton nombre exact de pixels mesuré entre Ordre 0 et H-alpha
    PIXELS_MESURES = 704.0 #278.0  
    
    # Paramètres de ta caméra
    PIXEL_NATIVE_UM = 4.8  # Taille native du pixel (µm)
    BINNING_UTILISE = 1    # Mettre 2 si ton image a été binnée en 2x2, sinon 1

    distance = compute_sa_distance_mm(
        num_pixels=PIXELS_MESURES,
        pixel_size_um=PIXEL_NATIVE_UM,
        binning=BINNING_UTILISE,
        lines_per_mm=200.0
    )
"""

# --- EXEMPLE D'UTILISATION ---
if __name__ == "__main__":
    # Nombre de pixels mesurés entre l'Ordre 0 et le creux d'absorption O2 à 7600 Å
    PIXELS_O2 = 820.0

    distance = compute_sa_distance_from_o2(
        num_pixels=PIXELS_O2,
        pixel_size_um=4.8,
        binning=1,       # 2 si binné 2x2 (soit 9.6 µm)
        lines_per_mm=200.0
    )

"""
