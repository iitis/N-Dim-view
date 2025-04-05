// CsvColumnAssignmentDialog.cpp
#include "CsvColumnAssignmentDialog.h"
#include <QLabel>
#include <QPushButton>
#include <QHBoxLayout>
#include <QMessageBox>

CsvColumnAssignmentDialog::CsvColumnAssignmentDialog(const QStringList& headers,
    const QVector<GroupDefinition>& userGroups,
    const QVector<QPair<int, int>>& defaults,
    QWidget* parent)
    : QDialog(parent), headers(headers) {
    setWindowTitle("Assign features to groups");
    mainLayout = new QVBoxLayout(this);
    formLayout = new QFormLayout;

    // Dodajemy domyślną grupę "pomiń" na początek
    groups.append(GroupDefinition( "<skip this feature>", {} ));
    groups += userGroups;

    buildUI();
    connectSignals();

    QPushButton* okButton = new QPushButton("OK");
    okButton->setDefault(true);
    connect(okButton, &QPushButton::clicked, this, [this]() {
        // Walidacja unikalności par grupa-etykieta
        int unnamed_cnt = 0;
        int spacial_cnt = 0;

        QSet<QString> used;
        for (const auto& widgets : rowWidgetsList) {
            int groupIdx = widgets.currentGroupIndex;
            int labelIdx = widgets.currentLabelIndex;

            if (groupIdx == 1) // note that there is <skip> group added as 0
                unnamed_cnt++;
            else if (groupIdx == 2)
                spacial_cnt++;

            if (groupIdx <= 0 || labelIdx < 0 || labelIdx >= groups[groupIdx].elementNames.size())
                continue; // pomiń lub brak etykiety

            QString key = QString::number(groupIdx) + ":" + groups[groupIdx].elementNames[labelIdx];
            if (used.contains(key)) {
                QMessageBox::warning(this, "N-Dim-view plugin error", "The same identifier has been assigned more than once in the same group.");
                return; // nie zamykamy dialogu
            }
            used.insert(key);
        }

        if ((unnamed_cnt < 1) && (spacial_cnt < 4)) {
            QMessageBox::warning(this, "N-Dim-view plugin error", "At least 1 element in the unnamed group or at least 4 elements in the spacial group must be assigned for the PCA to function properly.");
            return; // nie zamykamy dialogu
        }

        accept();
        });

    mainLayout->addLayout(formLayout);
    mainLayout->addWidget(okButton);

    auto defs = defaults;

    if (defs.empty()) {
        for (int i = 0; i < rowWidgetsList.size(); i++) {
            defs.push_back(QPair<int, int>(0, -1));
        }
    }

    int size = std::min(rowWidgetsList.size(), defs.size());

    for (int i = 0; i < size; i++) {
        int group = defs[i].first;
        rowWidgetsList[i].groupCombo->setCurrentIndex(group+1);
        if (group >= 0) {
            int label = defs[i].second;
            if (!userGroups[group].elementNames.empty() && label >= 0) {
                rowWidgetsList[i].labelCombo->setCurrentIndex(label);
            }
        }
    }
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
        //labelLayout->addWidget(new QLabel("Label:"));
        labelLayout->addWidget(labelCombo);
        labelLayout->setContentsMargins(0, 0, 0, 0);
        labelWidget->setVisible(false);

        QWidget* container = new QWidget;
        QHBoxLayout* vbox = new QHBoxLayout(container);
        vbox->addWidget(groupCombo);
        vbox->addWidget(labelWidget);
        vbox->setContentsMargins(0, 0, 0, 0);

        formLayout->addRow(headers[row], container);
        rowWidgetsList.append({ groupCombo, labelCombo, labelWidget });

        // Domyślnie ustawiamy grupę "<pomiń>"
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

        assignment.groupIndex = (widgets.currentGroupIndex < 0) ? -1 : widgets.currentGroupIndex - 1; // pomiń = -1, reszta od 0

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
