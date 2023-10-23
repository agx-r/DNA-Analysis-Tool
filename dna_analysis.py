import sys
from PyQt5.QtWidgets import QApplication, QWidget, QLabel, QLineEdit, QPushButton, QVBoxLayout, QMessageBox
from bioinformatics import calculate_gc_content, calculate_molecular_weight, reverse_complement

class DNAAnalysisWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("DNA Analysis Tool")
        self.initUI()

    def initUI(self):
        # Labels
        sequence_label = QLabel("DNA Sequence:")
        gc_label = QLabel("GC Content:")
        mw_label = QLabel("Molecular Weight:")
        rc_label = QLabel("Reverse Complement:")

        # Text fields
        self.sequence_field = QLineEdit()
        self.sequence_field.setPlaceholderText("Enter DNA sequence")
        self.sequence_field.setMaxLength(50)

        self.gc_content_field = QLineEdit()
        self.gc_content_field.setReadOnly(True)

        self.molecular_weight_field = QLineEdit()
        self.molecular_weight_field.setReadOnly(True)

        self.reverse_complement_field = QLineEdit()
        self.reverse_complement_field.setReadOnly(True)

        # Buttons
        analyze_button = QPushButton("Analyze")
        analyze_button.clicked.connect(self.analyze_sequence)

        # Layout
        layout = QVBoxLayout()
        layout.addWidget(sequence_label)
        layout.addWidget(self.sequence_field)
        layout.addWidget(gc_label)
        layout.addWidget(self.gc_content_field)
        layout.addWidget(mw_label)
        layout.addWidget(self.molecular_weight_field)
        layout.addWidget(rc_label)
        layout.addWidget(self.reverse_complement_field)
        layout.addWidget(analyze_button)

        self.setLayout(layout)

    def analyze_sequence(self):
        sequence = self.sequence_field.text().strip().upper()

        if not sequence:
            QMessageBox.warning(self, "Input Error", "Please enter a DNA sequence.")
            return

        try:
            gc_content = calculate_gc_content(sequence)
            molecular_weight = calculate_molecular_weight(sequence)
            reverse_comp = reverse_complement(sequence)

            self.gc_content_field.setText(f"{gc_content:.2f}%")
            self.molecular_weight_field.setText(f"{molecular_weight:.2f}")
            self.reverse_complement_field.setText(reverse_comp)
        except Exception as e:
            QMessageBox.critical(self, "Error", str(e))

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = DNAAnalysisWindow()
    window.show()
    sys.exit(app.exec_())
