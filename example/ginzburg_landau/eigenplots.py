import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from scipy.linalg import expm, svdvals


def plot_modal_analysis(savename):
    # > Load eigenspectrum and eigenvectors.
    eigvals = np.load("eigenspectrum.npy")
    direct_mode = np.load("eigenvectors.npy")

    # > Mesh.
    nx = len(direct_mode)
    x = np.linspace(-100, 100, nx+2)[1:-1]

    # > Analytic eigenspectrum.
    nu = 2 + 1j*0.2
    gamma = 1 - 1j
    mu_0 = 0.38
    c_mu = 0.2
    mu_2 = -0.01
    h = np.sqrt(-2*mu_2*gamma)

    def true_spec(n): return mu_0 - c_mu**2 - nu**2/(4*gamma) - (n + 0.5)*h

    # --------------------------
    # -----     FIGURE     -----
    # --------------------------

    fig, axes = plt.subplot_mosaic(
        [["spectrum", "mode", "mode"]], figsize=(9, 3), layout="constrained")

    # True eigenspectrum.
    axes["spectrum"].plot(
        true_spec(np.arange(40)).imag, true_spec(np.arange(40)).real,
        "o", color="black", mfc="none",
        label="Analytic"
    )

    # Direct eigenspectrum.
    axes["spectrum"].plot(
        eigvals[:, 1], eigvals[:, 0],
        ".", color="dodgerblue", ms=5,
        label="LightKrylov"
    )

    # Aesthetics.
    axes["spectrum"].legend(
        loc="lower center", bbox_to_anchor=(0.5, 1.05), ncols=2)
    axes["spectrum"].set(xlabel=r"$\omega$", xlim=(-2, 2))
    axes["spectrum"].set(ylabel=r"$\sigma$", ylim=(-4, 1))

    rect = Rectangle(
        (-2, 0),  # Upper left corner.
        4, -4,   # Width / Height
        edgecolor="black", lw=0.5,
        facecolor="lightgray", alpha=0.5,
    )

    axes["spectrum"].add_patch(rect)

    axes["mode"].plot(x, direct_mode[:nx, 0].real, lw=1,
                      c="dodgerblue", label="Leading eigenvector")
    axes["mode"].plot(x, direct_mode[:nx, 0].imag,
                      lw=1, c="dodgerblue", ls="--")
    axes["mode"].plot(x, abs(direct_mode[:nx, 0]), lw=0.5, c="k")
    axes["mode"].plot(x, -abs(direct_mode[:nx, 0]), lw=0.5, c="k")

    axes["mode"].set(yticks=[])
    axes["mode"].set(xlim=(-40, 40), xlabel=r"$x$")

    axes["mode"].legend(loc="lower center",
                        bbox_to_anchor=(0.5, 1.05), ncols=2)

    plt.savefig(savename, bbox_inches="tight", dpi=1200)


if __name__ == "__main__":
    plot_modal_analysis("ginzburg_landau_modal_analysis.png")
    plt.show()
