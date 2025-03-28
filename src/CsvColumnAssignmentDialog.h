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
};

struct ColumnAssignment {
    int featureIndex; // indeks cechy - wymiaru w pliku Ÿród³owym
    QString featureName; // nazwa cechy - jak wpisano w pierwszej linii pliku csv
    int groupIndex; // indeks grupy - lub -1 jesli cecha nie ma byæ dalej analizowana
    QString groupName; // nazwa tej grupy - taka jak sobie zdefiniujesz lub pusty QString jeœli by³o -1
    std::optional<int> label_id; // indeks wybranej etykiety dla cechy w obrêbie grupy (jesli grupa ma etykiety) 
    std::optional<QString> label_name; // nazwa tej etykiety taka jak sopbie zdefiniujesz (o ile oczywiœcie grupa ma etykiety)
};

class CsvColumnAssignmentDialog : public QDialog {
    Q_OBJECT

public:
    CsvColumnAssignmentDialog(const QStringList& headers,
        const QVector<GroupDefinition>& groups,
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
