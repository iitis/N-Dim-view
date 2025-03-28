// CsvColumnAssignmentDialog.cpp
#include "CsvColumnAssignmentDialog.h"
#include <QLabel>
#include <QPushButton>
#include <QHBoxLayout>
#include <QMessageBox>

CsvColumnAssignmentDialog::CsvColumnAssignmentDialog(const QStringList& headers,
    const QVector<GroupDefinition>& userGroups,
    QWidget* parent)
    : QDialog(parent), headers(headers) {
    setWindowTitle("Assign features to groups");
    mainLayout = new QVBoxLayout(this);
    formLayout = new QFormLayout;

    // Dodajemy domy�ln� grup� "pomi�" na pocz�tek
    groups.append({ "<skip this feature>", {} });
    groups += userGroups;

    buildUI();
    connectSignals();

    QPushButton* okButton = new QPushButton("OK");
    okButton->setDefault(true);
    connect(okButton, &QPushButton::clicked, this, [this]() {
        // Walidacja unikalno�ci par grupa-etykieta
        QSet<QString> used;
        for (const auto& widgets : rowWidgetsList) {
            int groupIdx = widgets.currentGroupIndex;
            int labelIdx = widgets.currentLabelIndex;

            if (groupIdx <= 0 || labelIdx < 0 || labelIdx >= groups[groupIdx].elementNames.size())
                continue; // pomi� lub brak etykiety

            QString key = QString::number(groupIdx) + ":" + groups[groupIdx].elementNames[labelIdx];
            if (used.contains(key)) {
                QMessageBox::warning(this, "Error", "The same identifier has been assigned more than once in the same group.");
                return; // nie zamykamy dialogu
            }
            used.insert(key);
        }
        accept();
        });

    mainLayout->addLayout(formLayout);
    mainLayout->addWidget(okButton);
}


void CsvColumnAssignmentDialog::buildUI() {
    for (int row = 0; row < headers.size(); ++row) {
        QComboBox* groupCombo = new QComboBox;
        for (const auto& group : groups)
            groupCombo->addItem(group.name);

        QComboBox* labelCombo = new QComboBox;
        labelCombo->addItem("<none>");
        labelCombo->setEnabled(false);

        QWidget* labelWidget = new QWidget;
        QHBoxLayout* labelLayout = new QHBoxLayout(labelWidget);
        labelLayout->addWidget(new QLabel("Label:"));
        labelLayout->addWidget(labelCombo);
        labelLayout->setContentsMargins(0, 0, 0, 0);
        labelWidget->setVisible(false);

        QWidget* container = new QWidget;
        QVBoxLayout* vbox = new QVBoxLayout(container);
        vbox->addWidget(groupCombo);
        vbox->addWidget(labelWidget);
        vbox->setContentsMargins(0, 0, 0, 0);

        formLayout->addRow(headers[row], container);
        rowWidgetsList.append({ groupCombo, labelCombo, labelWidget });

        // Domy�lnie ustawiamy grup� "<pomi�>"
        groupCombo->setCurrentIndex(0);
    }
}

void CsvColumnAssignmentDialog::updateLabelCombo(RowWidgets& widgets, const GroupDefinition& group) {
    widgets.labelCombo->blockSignals(true);
    widgets.labelCombo->clear();

    if (!group.elementNames.isEmpty()) {
        widgets.labelCombo->addItems(group.elementNames);
        widgets.labelCombo->setEnabled(true);
        widgets.labelWidget->setVisible(true);
        widgets.labelCombo->setCurrentIndex(0);
        widgets.currentLabelIndex = 0;
    }
    else {
        widgets.labelCombo->addItem("<none>");
        widgets.labelCombo->setEnabled(false);
        widgets.labelWidget->setVisible(false);
        widgets.currentLabelIndex = -1;
    }

    widgets.labelCombo->blockSignals(false);
}

void CsvColumnAssignmentDialog::connectSignals() {
    for (int i = 0; i < rowWidgetsList.size(); ++i) {
        auto& row = rowWidgetsList[i];

        connect(row.groupCombo, QOverload<int>::of(&QComboBox::currentIndexChanged),
            [this, i](int newGroupIndex) {
                auto& widgets = rowWidgetsList[i];
                widgets.currentGroupIndex = newGroupIndex;
                updateLabelCombo(widgets, groups[newGroupIndex]);
            });

        connect(row.labelCombo, QOverload<int>::of(&QComboBox::currentIndexChanged),
            [this, i](int newLabelIndex) {
                auto& widgets = rowWidgetsList[i];
                widgets.currentLabelIndex = newLabelIndex;
            });
    }
}

QVector<ColumnAssignment> CsvColumnAssignmentDialog::getAssignments() const {
    QVector<ColumnAssignment> result;
    for (int i = 0; i < rowWidgetsList.size(); ++i) {
        const auto& widgets = rowWidgetsList[i];
        ColumnAssignment assignment;
        assignment.featureIndex = i;
        assignment.featureName = headers[i];

        assignment.groupIndex = (widgets.currentGroupIndex < 0) ? -1 : widgets.currentGroupIndex - 1; // pomi� = -1, reszta od 0

        assignment.groupName = (assignment.groupIndex < 0) ? QString() : groups[widgets.currentGroupIndex].name;

        if (widgets.currentGroupIndex > 0 && !groups[widgets.currentGroupIndex].elementNames.isEmpty()) {
            assignment.label_id = widgets.currentLabelIndex;
            if (widgets.currentLabelIndex >= 0 && widgets.currentLabelIndex < groups[widgets.currentGroupIndex].elementNames.size())
                assignment.label_name = groups[widgets.currentGroupIndex].elementNames[widgets.currentLabelIndex];
        }

        result.append(assignment);
    }
    return result;
}
