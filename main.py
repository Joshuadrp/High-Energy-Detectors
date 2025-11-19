# main.py
import subprocess
import sys

def main():
    print("\nGAMMA-RAY DETECTOR ANALYSIS")
    print("\nSelect a detector to analyze:")
    print("\n  1. NAITI(Thallium-doped Sodium Iodide)")
    print("  2. BGO (Bismuth Germanate)")
    print("  3. CdTe (Cadmium Telluride)")
    print("  4. Exit")

    choice = input("\nEnter your choice (1-4): ").strip()

    if choice == '1':
        print("\nRunning NAITI analysis...\n")
        subprocess.run([sys.executable, "NAITI.py"])

    elif choice == '2':
        print("\nRunning BGO analysis...\n")
        subprocess.run([sys.executable, "BGO.py"])

    elif choice == '3':
        print("\nRunning CdTe analysis...\n")
        subprocess.run([sys.executable, "CDTE.py"])

    elif choice == '4':
        print("\nGoodbye!\n")
        sys.exit(0)

    else:
        print("\nInvalid choice! Please enter 1-4.")
        main()

if __name__ == "__main__":
    main()