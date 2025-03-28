#pragma once

#include <QVector>
#include <Eigen/Dense>

#include <QMap>
#include <QInputDialog>
#include <QDialog>
#include <QRegularExpression>

#include <optional>


class CsvReader
{
public:
    struct CsvData {
        QStringList headers;
        QVector<Eigen::VectorXd> samples;
    };

    static CsvData loadFile(const QString & = QString());
    static CsvData loadCsvToEigenVectors(const QString& filePath, QRegularExpression regsep = QRegularExpression("[,;\t]"));
    static Eigen::MatrixXd convertCsvDataToMatrix(const CsvData& csv_data);
};


