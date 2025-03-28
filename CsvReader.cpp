#include "CsvReader.h"

#include <QFileDialog>
#include <QDir>
#include <QFile>
#include <QTextStream>

#include <QStringList>
#include <QDebug>


CsvReader::CsvData CsvReader::loadFile(const QString& path) {
    QString fileName = path;
    if (fileName.isEmpty() || !QFileInfo(fileName).exists()) {
        qInfo() << "path is empty. Let's to use the FileDialog now.";

        fileName = QDir::toNativeSeparators(
            QFileDialog::getOpenFileName(0, QString::fromUtf8("Select .csv/.dat file"), "c:/K3/Wielowymiar/VinhoVerde/*.*", "ALL (*.*);;CSV (*.csv);;DAT (*.dat)")
        );

        if (fileName.isEmpty() || !QFileInfo(fileName).exists()) {
            return CsvData();
        }
    }

    return loadCsvToEigenVectors(fileName);
}


CsvReader::CsvData CsvReader::loadCsvToEigenVectors(const QString& filePath, QRegularExpression regsep) {
    CsvData result;

    QFile file(filePath);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        qWarning() << "Can't open file:" << filePath;
        return result;
    }

    QTextStream in(&file);
    bool firstLine = true;
    int dimension = 0;

    while (!in.atEnd()) {
        QString line = in.readLine().trimmed();
        if (line.isEmpty())
            continue;

        QStringList parts = line.split(regsep);

        if (firstLine) {
            bool hasText = false;
            for (QString part : parts) {
                part.remove(QRegularExpression("^\"|\"$")); // usuń cudzysłowy
                bool ok = false;
                part.toDouble(&ok);
                if (!ok) {
                    hasText = true;
                    break;
                }
            }

            dimension = parts.size();
            firstLine = false;

            if (hasText) {
                // Mamy teksty – to są nagłówki
                for (QString& header : parts)
                    header = header.remove(QRegularExpression("^\"|\"$"));
                result.headers = parts;
                continue;
            }
            else {
                // Brak nagłówków – generujemy własne
                result.headers.clear();
                for (int i = 0; i < dimension; ++i)
                    result.headers << QString("Feature %1").arg(i + 1);
                // spadamy dalej i przetwarzamy tę linię jako dane
            }
        }

        if (parts.size() != dimension) {
            qWarning() << "The number of columns does not match the inline header:" << line;
            continue;
        }

        Eigen::VectorXd vec(dimension);
        for (int i = 0; i < dimension; ++i) {
            bool ok = false;
            double value = parts[i].toDouble(&ok);
            if (!ok) {
                qWarning() << "Incorrect value in line:" << line;
                vec[i] = 0.0;
            }
            else {
                vec[i] = value;
            }
        }
        result.samples.append(vec);
    }

    return result;
}

Eigen::MatrixXd CsvReader::convertCsvDataToMatrix(const CsvData& csv_data)
{
    if (csv_data.samples.isEmpty()) return Eigen::MatrixXd();

    int numRows = csv_data.samples.size();         // liczba próbek (wiersze)
    int numCols = csv_data.samples[0].size();      // liczba wymiarów (kolumny)

    Eigen::MatrixXd mat(numCols, numRows); // odwrotnie niż wcześniej

    for (int row = 0; row < numRows; ++row) {
        for (int col = 0; col < numCols; ++col) {
            mat(col, row) = csv_data.samples[row](col);
        }
    }

    return mat;
}


