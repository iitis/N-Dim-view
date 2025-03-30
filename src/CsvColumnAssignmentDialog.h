// CsvColumnAssignmentDialog.h
#pragma once

#include <QDialog>
#include <QComboBox>
#include <QFormLayout>
#include <QVector>
#include <QStringList>
#include <QVBoxLayout>
#include <QMap>
#include <QSet>
#include <optional>

struct GroupDefinition {
    QString name;
    QStringList elementNames; // pusty = brak ograniczenia

    GroupDefinition(const QString& n, QStringList c = {})
    : name(n), elementNames(std::move(c)) {}
};

struct ColumnAssignment {
    int featureIndex; // indeks cechy - wymiaru w pliku źródłowym
    QString featureName; // nazwa cechy - jak wpisano w pierwszej linii pliku csv
    int groupIndex; // indeks grupy - lub -1 jesli cecha nie ma być dalej analizowana
    QString groupName; // nazwa tej grupy - taka jak sobie zdefiniujesz lub pusty QString jeśli było -1
    std::optional<int> label_id; // indeks wybranej etykiety dla cechy w obrębie grupy (jesli grupa ma etykiety) 
    std::optional<QString> label_name; // nazwa tej etykiety taka jak sopbie zdefiniujesz (o ile oczywiście grupa ma etykiety)
};

class CsvColumnAssignmentDialog : public QDialog {
    Q_OBJECT

public:
    CsvColumnAssignmentDialog(const QStringList& headers,
        const QVector<GroupDefinition>& groups,
        const QVector<QPair<int,int>>& defaults = QVector<QPair<int, int>>(),
        QWidget* parent = nullptr);

    QVector<ColumnAssignment> getAssignments() const;

private:
    QStringList headers;
    QVector<GroupDefinition> groups;
    QFormLayout* formLayout;
    QVBoxLayout* mainLayout;

    struct RowWidgets {
        QComboBox* groupCombo;
        QComboBox* labelCombo;
        QWidget* labelWidget;
        int currentGroupIndex = -1;
        int currentLabelIndex = -1;
    };

    QVector<RowWidgets> rowWidgetsList;

    void buildUI();
    void connectSignals();
    void updateLabelCombo(struct RowWidgets& widgets, const GroupDefinition& group);

};
