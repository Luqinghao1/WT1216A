#ifndef FITTINGWIDGET_H
#define FITTINGWIDGET_H

#include <QWidget>
#include <QFutureWatcher>
#include <QApplication>
#include <QMouseEvent>
#include <QDialog>
#include <QFormLayout>
#include <QCheckBox>
#include <QLineEdit>
#include <QSpinBox>
#include <QDialogButtonBox>
#include <QPushButton>
#include <QTableWidget>
#include <QComboBox>
#include <QLabel>
#include <QGroupBox>
#include <QProgressBar>
#include <QSplitter>
#include <QHeaderView>

// 引入项目通用模块
#include "modelmanager.h"
#include "mousezoom.h"
#include "chartsetting1.h"
#include "modelselect.h" // 引入模型选择对话框

namespace Ui {
class FittingWidget;
}

// 拟合参数结构体
struct FitParameter {
    QString name;
    QString displayName;
    QString symbol;
    double value;
    double min;
    double max;
    bool isFit;
    QString unit;
};

class FittingWidget : public QWidget
{
    Q_OBJECT

public:
    explicit FittingWidget(QWidget *parent = nullptr);
    ~FittingWidget();

    // 设置模型管理器
    void setModelManager(ModelManager* m);
    // 设置实测数据
    void setObservedData(const QVector<double>& t, const QVector<double>& p, const QVector<double>& d);

signals:
    // 迭代更新信号（用于更新图表和误差显示）
    void sigIterationUpdated(double error, const QMap<QString,double>& currentParams,
                             const QVector<double>& t, const QVector<double>& p, const QVector<double>& d);
    // 进度信号
    void sigProgress(int percent);
    // 拟合完成信号
    void fittingCompleted(ModelManager::ModelType type, const QMap<QString,double>& finalParams);

private slots:
    // --- 界面按钮响应槽 ---
    void on_btnLoadData_clicked();
    void on_btnRunFit_clicked();
    void on_btnStop_clicked();
    void on_btnResetParams_clicked(); // 重置参数（根据当前模型类型）
    void on_btnImportModel_clicked();
    void on_btnResetView_clicked();
    void on_btnExportData_clicked();
    void on_btnExportChart_clicked();
    void on_btnChartSettings_clicked();

    // 修改：点击模型选择按钮的槽函数
    void on_btn_modelSelect_clicked();

    // --- 逻辑处理槽 ---
    void onIterationUpdate(double err, const QMap<QString,double>& p,
                           const QVector<double>& t, const QVector<double>& p_curve, const QVector<double>& d_curve);
    void onFitFinished();

private:
    Ui::FittingWidget *ui;
    ModelManager* m_modelManager;
    MouseZoom* m_plot;
    QCPTextElement* m_plotTitle;

    bool m_isFitting;
    bool m_stopRequested;

    // 当前选中的模型类型
    ModelManager::ModelType m_currentModelType;

    // 数据与参数
    QVector<double> m_obsTime;
    QVector<double> m_obsPressure;
    QVector<double> m_obsDerivative;
    QList<FitParameter> m_parameters;
    QFutureWatcher<void> m_watcher;

    // --- 内部辅助函数 ---
    void setupPlot();
    void initializeDefaultModel(); // 初始化默认模型状态

    void loadParamsToTable();
    void updateParamsFromTable();
    void updateModelCurve();

    void getParamDisplayInfo(const QString& key, QString& outName, QString& outSymbol, QString& outUnicodeSymbol, QString& outUnit);
    QStringList getParamOrder(ModelManager::ModelType type);
    QStringList parseLine(const QString& line);
    void plotCurves(const QVector<double>& t, const QVector<double>& p, const QVector<double>& d, bool isModel);

    // --- 拟合算法相关 ---
    void runOptimizationTask(ModelManager::ModelType modelType, QList<FitParameter> fitParams, double weight);
    void runLevenbergMarquardtOptimization(ModelManager::ModelType modelType, QList<FitParameter> params, double weight);
    QVector<double> calculateResiduals(const QMap<QString,double>& params, ModelManager::ModelType modelType, double weight);
    double calculateSumSquaredError(const QVector<double>& residuals);
    QVector<QVector<double>> computeJacobian(const QMap<QString,double>& params,
                                             const QVector<double>& baseResiduals,
                                             const QVector<int>& fitIndices,
                                             ModelManager::ModelType modelType,
                                             const QList<FitParameter>& currentFitParams,
                                             double weight);
    QVector<double> solveLinearSystem(const QVector<QVector<double>>& A, const QVector<double>& b);
};

// 数据加载对话框
class FittingDataLoadDialog : public QDialog {
    Q_OBJECT
public:
    FittingDataLoadDialog(const QList<QStringList>& previewData, QWidget* parent=nullptr);
    int getTimeColumnIndex() const;
    int getPressureColumnIndex() const;
    int getDerivativeColumnIndex() const;
    int getSkipRows() const;
    int getPressureDataType() const;
private slots:
    void validateSelection();
private:
    QTableWidget* m_previewTable;
    QComboBox* m_comboTime;
    QComboBox* m_comboPressure;
    QComboBox* m_comboDeriv;
    QComboBox* m_comboSkipRows;
    QComboBox* m_comboPressureType;
};

#endif // FITTINGWIDGET_H
