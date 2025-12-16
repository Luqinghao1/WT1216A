#include "modelmanager.h"
#include "modelselect.h"
#include "modelwidget1.h"
#include "modelwidget2.h"
#include "modelwidget3.h"
#include "modelwidget4.h"
#include "modelwidget5.h"
#include "modelwidget6.h"

#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QGroupBox>
#include <QDebug>
#include <cmath>

ModelManager::ModelManager(QWidget* parent)
    : QObject(parent), m_mainWidget(nullptr), m_btnSelectModel(nullptr), m_modelStack(nullptr)
    , m_modelWidget1(nullptr), m_modelWidget2(nullptr), m_modelWidget3(nullptr)
    , m_modelWidget4(nullptr), m_modelWidget5(nullptr), m_modelWidget6(nullptr)
    , m_currentModelType(Model_1)
{
}

ModelManager::~ModelManager() {}

void ModelManager::initializeModels(QWidget* parentWidget)
{
    if (!parentWidget) return;
    createMainWidget();
    setupModelSelection();

    m_modelStack = new QStackedWidget(m_mainWidget);

    m_modelWidget1 = new ModelWidget1(m_modelStack);
    m_modelWidget2 = new ModelWidget2(m_modelStack);
    m_modelWidget3 = new ModelWidget3(m_modelStack);
    m_modelWidget4 = new ModelWidget4(m_modelStack);
    m_modelWidget5 = new ModelWidget5(m_modelStack);
    m_modelWidget6 = new ModelWidget6(m_modelStack);

    m_modelStack->addWidget(m_modelWidget1); // Index 0
    m_modelStack->addWidget(m_modelWidget2); // Index 1
    m_modelStack->addWidget(m_modelWidget3); // Index 2
    m_modelStack->addWidget(m_modelWidget4); // Index 3
    m_modelStack->addWidget(m_modelWidget5); // Index 4
    m_modelStack->addWidget(m_modelWidget6); // Index 5

    m_mainWidget->layout()->addWidget(m_modelStack);
    connectModelSignals();

    switchToModel(Model_1);

    if (parentWidget->layout()) parentWidget->layout()->addWidget(m_mainWidget);
    else {
        QVBoxLayout* layout = new QVBoxLayout(parentWidget);
        layout->addWidget(m_mainWidget);
        parentWidget->setLayout(layout);
    }
}

void ModelManager::createMainWidget()
{
    m_mainWidget = new QWidget();
    QVBoxLayout* mainLayout = new QVBoxLayout(m_mainWidget);
    mainLayout->setContentsMargins(10, 5, 10, 10);
    mainLayout->setSpacing(5);
    m_mainWidget->setLayout(mainLayout);
}

void ModelManager::setupModelSelection()
{
    if (!m_mainWidget) return;
    QGroupBox* selectionGroup = new QGroupBox("模型选择", m_mainWidget);
    selectionGroup->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed);
    QHBoxLayout* selectionLayout = new QHBoxLayout(selectionGroup);
    selectionLayout->setContentsMargins(9, 9, 9, 9);

    QLabel* infoLabel = new QLabel("当前模型:", selectionGroup);

    m_btnSelectModel = new QPushButton("点击选择模型...", selectionGroup);
    m_btnSelectModel->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
    m_btnSelectModel->setMinimumHeight(30);
    m_btnSelectModel->setStyleSheet("text-align: left; padding-left: 10px; font-weight: bold;");

    connect(m_btnSelectModel, &QPushButton::clicked, this, &ModelManager::onSelectModelClicked);

    selectionLayout->addWidget(infoLabel);
    selectionLayout->addWidget(m_btnSelectModel);

    QVBoxLayout* mainLayout = qobject_cast<QVBoxLayout*>(m_mainWidget->layout());
    if (mainLayout) {
        mainLayout->addWidget(selectionGroup);
    }
}

void ModelManager::connectModelSignals()
{
    if (m_modelWidget1) connect(m_modelWidget1, &ModelWidget1::calculationCompleted, this, &ModelManager::onWidgetCalculationCompleted);
    if (m_modelWidget2) connect(m_modelWidget2, &ModelWidget2::calculationCompleted, this, &ModelManager::onWidgetCalculationCompleted);
}

void ModelManager::switchToModel(ModelType modelType)
{
    if (!m_modelStack) return;
    ModelType old = m_currentModelType;
    m_currentModelType = modelType;
    m_modelStack->setCurrentIndex((int)modelType);

    QString name = getModelTypeName(modelType);
    if (m_btnSelectModel) m_btnSelectModel->setText(name);

    emit modelSwitched(modelType, old);
}

void ModelManager::onSelectModelClicked()
{
    ModelSelect dlg(m_mainWidget);
    if (dlg.exec() == QDialog::Accepted) {
        QString code = dlg.getSelectedModelCode();
        if (code == "modelwidget1") switchToModel(Model_1);
        else if (code == "modelwidget2") switchToModel(Model_2);
        else if (code == "modelwidget3") switchToModel(Model_3);
        else if (code == "modelwidget4") switchToModel(Model_4);
        else if (code == "modelwidget5") switchToModel(Model_5);
        else if (code == "modelwidget6") switchToModel(Model_6);
        else { }
    }
}

QString ModelManager::getModelTypeName(ModelType type)
{
    switch (type) {
    case Model_1: return "压裂水平井复合页岩油模型1";
    case Model_2: return "压裂水平井复合页岩油模型2";
    case Model_3: return "压裂水平井复合页岩油模型3";
    case Model_4: return "压裂水平井复合页岩油模型4";
    case Model_5: return "压裂水平井复合页岩油模型5";
    case Model_6: return "压裂水平井复合页岩油模型6";
    default: return "未知模型";
    }
}

void ModelManager::onWidgetCalculationCompleted(const QString &t, const QMap<QString, double> &r) {
    emit calculationCompleted(t, r);
}

void ModelManager::setHighPrecision(bool high) {
    if (m_modelWidget1) m_modelWidget1->setHighPrecision(high);
    if (m_modelWidget2) m_modelWidget2->setHighPrecision(high);
}

QMap<QString, double> ModelManager::getDefaultParameters(ModelType type)
{
    QMap<QString, double> p;
    p.insert("cD", 0.001);
    p.insert("S", 0.01);
    p.insert("gamaD", 0.02);

    if (type == Model_1 || type == Model_2) {
        p.insert("kf", 1e-3);
        p.insert("km", 1e-4);
        p.insert("L", 1000.0);
        p.insert("Lf", 100.0);
        p.insert("LfD", 0.1);
        p.insert("rmD", 4.0);
        p.insert("omega1", 0.4);
        p.insert("omega2", 0.08);
        p.insert("lambda1", 1e-3);
        p.insert("gamaD", 0.02);

        if (type == Model_1) {
            // --- 修改：Model 1 默认参数 ---
            p.insert("cD", 0.01);
            p.insert("S", 1.0);
        } else {
            // Model 2
            p.insert("cD", 0.0);
            p.insert("S", 0.0);
        }
    }
    return p;
}

ModelCurveData ModelManager::calculateTheoreticalCurve(ModelType type, const QMap<QString, double>& params, const QVector<double>& providedTime)
{
    if (type == Model_1 && m_modelWidget1) {
        return m_modelWidget1->calculateTheoreticalCurve(params, providedTime);
    }
    else if (type == Model_2 && m_modelWidget2) {
        return m_modelWidget2->calculateTheoreticalCurve(params, providedTime);
    }
    return ModelCurveData();
}

QVector<double> ModelManager::generateLogTimeSteps(int count, double startExp, double endExp) {
    QVector<double> t;
    t.reserve(count);
    for (int i = 0; i < count; ++i) {
        double exponent = startExp + (endExp - startExp) * i / (count - 1);
        t.append(pow(10.0, exponent));
    }
    return t;
}

void ModelManager::setObservedData(const QVector<double>& t, const QVector<double>& p, const QVector<double>& d)
{
    m_cachedObsTime = t;
    m_cachedObsPressure = p;
    m_cachedObsDerivative = d;
}

void ModelManager::getObservedData(QVector<double>& t, QVector<double>& p, QVector<double>& d) const
{
    t = m_cachedObsTime;
    p = m_cachedObsPressure;
    d = m_cachedObsDerivative;
}

bool ModelManager::hasObservedData() const
{
    return !m_cachedObsTime.isEmpty();
}
