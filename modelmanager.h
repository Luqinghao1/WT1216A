#ifndef MODELMANAGER_H
#define MODELMANAGER_H

#include <QObject>
#include <QWidget>
#include <QStackedWidget>
#include <QPushButton>
#include <QMap>
#include <QVector>
#include <tuple>

// 定义模型曲线数据类型: <时间(t), 压力(pD), 导数(dpD)>
typedef std::tuple<QVector<double>, QVector<double>, QVector<double>> ModelCurveData;

// 前置声明
class ModelWidget1;
class ModelWidget2;
class ModelWidget3;
class ModelWidget4;
class ModelWidget5;
class ModelWidget6;

class ModelManager : public QObject
{
    Q_OBJECT

public:
    // 6 个模型
    enum ModelType {
        Model_1 = 0, // 压裂水平井复合页岩油模型1 (变井储+无穷大)
        Model_2 = 1, // 压裂水平井复合页岩油模型2 (恒定井储+无穷大)
        Model_3 = 2, // 压裂水平井复合页岩油模型3 (变井储+封闭)
        Model_4 = 3, // 压裂水平井复合页岩油模型4 (恒定井储+封闭)
        Model_5 = 4, // 压裂水平井复合页岩油模型5 (变井储+定压)
        Model_6 = 5  // 压裂水平井复合页岩油模型6 (恒定井储+定压)
    };
    Q_ENUM(ModelType)

    explicit ModelManager(QWidget* parent = nullptr);
    ~ModelManager();

    void initializeModels(QWidget* parentWidget);
    QWidget* getMainWidget() const { return m_mainWidget; }
    void switchToModel(ModelType modelType);
    ModelType getCurrentModelType() const { return m_currentModelType; }

    static QString getModelTypeName(ModelType type);

    // 设置高精度模式 (用于拟合时降低精度加速)
    void setHighPrecision(bool high);

    // 获取指定模型的默认参数
    QMap<QString, double> getDefaultParameters(ModelType type);

    // 统一计算接口，内部代理给具体的 ModelWidget
    ModelCurveData calculateTheoreticalCurve(ModelType type,
                                             const QMap<QString, double>& params,
                                             const QVector<double>& providedTime = QVector<double>());

    // --- 修复：加回生成对数时间步长的静态辅助函数 ---
    static QVector<double> generateLogTimeSteps(int count, double startExp, double endExp);

    // 实测数据持久化接口
    void setObservedData(const QVector<double>& t, const QVector<double>& p, const QVector<double>& d);
    void getObservedData(QVector<double>& t, QVector<double>& p, QVector<double>& d) const;
    bool hasObservedData() const;

signals:
    void modelSwitched(ModelType newType, ModelType oldType);
    void calculationCompleted(const QString& title, const QMap<QString, double>& results);

private slots:
    // 弹出模型选择对话框
    void onSelectModelClicked();
    // 接收子 Widget 计算完成信号并转发
    void onWidgetCalculationCompleted(const QString& t, const QMap<QString, double>& r);

private:
    void createMainWidget();
    void setupModelSelection(); // 初始化选择按钮
    void connectModelSignals();

    QWidget* m_mainWidget;
    QPushButton* m_btnSelectModel;
    QStackedWidget* m_modelStack;

    // 六个模型 Widget 指针
    ModelWidget1* m_modelWidget1;
    ModelWidget2* m_modelWidget2;
    ModelWidget3* m_modelWidget3;
    ModelWidget4* m_modelWidget4;
    ModelWidget5* m_modelWidget5;
    ModelWidget6* m_modelWidget6;

    ModelType m_currentModelType;

    // 缓存的实测数据
    QVector<double> m_cachedObsTime;
    QVector<double> m_cachedObsPressure;
    QVector<double> m_cachedObsDerivative;
};

#endif // MODELMANAGER_H
