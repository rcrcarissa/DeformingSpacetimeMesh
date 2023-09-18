import pickle
import matplotlib.pyplot as plt

num_upsampled = 1024

if __name__ == "__main__":
    # with open("psnr_straight.pickle", "rb") as f:
    #     psnr_straight = pickle.load(f)
    # with open("psnr_deformed.pickle", "rb") as f:
    #     psnr_deformed = pickle.load(f)
    with open("psnr_straight_denoise.pickle", "rb") as f:
        psnr_straight_denoise = pickle.load(f)
    with open("psnr_deformed_denoise.pickle", "rb") as f:
        psnr_deformed_denoise = pickle.load(f)
    x = range(num_upsampled)

    plt.rcParams['font.sans-serif'] = ['Arial']  # font
    plt.rcParams['axes.unicode_minus'] = False  # show minus sign

    plt.figure(figsize=(10, 6.5))
    plt.grid(linestyle="--")  # dashed line as grid
    ax = plt.gca()
    ax.spines['top'].set_visible(False)  # invisible top border
    ax.spines['right'].set_visible(False)  # invisible right border

    plt.plot(x, psnr_straight_denoise, color="blue", label="SL", linewidth=1.5)
    plt.plot(x, psnr_deformed_denoise, color="green", label="PL", linewidth=1.5)

    group_labels = ['0', r'$\frac{1}{2}\pi$', r'$\pi$', r'$\frac{3}{2}\pi$', r'$2\pi$']
    plt.xticks([0, num_upsampled / 4, num_upsampled / 2, num_upsampled / 4 * 3, num_upsampled], group_labels,
               fontsize=20, fontweight='bold')
    plt.yticks(fontsize=20, fontweight='bold')
    # plt.title("example", fontsize=12, fontweight='bold')
    plt.xlabel("$\phi$", fontsize=25, fontweight='bold', loc='right')
    plt.ylabel("PSNR", fontsize=25, fontweight='bold', loc='top')
    # plt.xlim(0.9, 6.1)
    plt.ylim(0, 60)

    plt.legend(loc=0, numpoints=1)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=20, fontweight='bold')

    # plt.show()
    plt.savefig("denoised_comp.svg")
