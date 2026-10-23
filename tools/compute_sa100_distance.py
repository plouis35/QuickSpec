import logging

# Constantes physiques
LAMBDA_O2_NM = 760.0       # 7600 Å (bande d'absorption atmosphérique de l'oxygène)
LAMBDA_HALPHA_NM = 656.28  # Longueur d'onde H-alpha (6562.8 Å)

def compute_sa_distance_mm(lambda_ref: float,
                           num_pixels: float, 
                           pixel_size_um: float = 4.8, 
                           binning: int = 1, 
                           lines_per_mm: float = 200.0) -> float:
    """
    Calcule la distance réseau SA -> capteur en mm à partir de la raie de réference.

    Args:
        ref_line (float): lambda (en nm) de la raie de référence utilisée
        num_pixels (float): Nombre de pixels entre l'Ordre 0 et la raie de référence
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
    lambda_halpha_mm = lambda_ref * 1e-6

    distance_mm = (d_mm * x_sensor_mm) / lambda_halpha_mm

    # Informations complémentaires
    dispersion_AA_per_px = (lambda_ref * 10.0) / num_pixels

    print(f"--- RÉSULTATS ---")
    print(f"Distance mesurée Ordre 0 -> raie   : {num_pixels:.2f} px")
    print(f"Taille de pixel effective          : {effective_pixel_um:.2f} µm (Binning {binning}x{binning})")
    print(f"Dispersion mesurée sur le spectre  : {dispersion_AA_per_px:.2f} Å/px")
    print(f"-----------------------------------")
    print(f"Distance SA -> Capteur             : {distance_mm:.2f} mm")

    return distance_mm

if __name__ == "__main__":
    PIXELS_MESURES = 704.0 #278.0  
    
    PIXEL_NATIVE_UM = 2.4  # Taille native du pixel de la caméra (µm)
    BINNING_UTILISE = 2    

    distance = compute_sa_distance_mm(
        lambda_ref=LAMBDA_HALPHA_NM,
        num_pixels=PIXELS_MESURES,
        pixel_size_um=PIXEL_NATIVE_UM,
        binning=BINNING_UTILISE,
        lines_per_mm=200.0
    )

