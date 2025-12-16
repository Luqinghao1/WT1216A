#ifndef MODELWIDGET2_H
#define MODELWIDGET2_H

#include <QWidget>
#include <QMap>
#include <QVector>
#include <functional>
#include <tuple>
#include <QLineEdit>

#include "mousezoom.h"
#include "chartsetting1.h"

typedef std::tuple<QVector<double>, QVector<double>, QVector<double>> ModelCurveData;

namespace Ui {
class ModelWidget2;
}

class ModelWidget2 : public QWidget
{
    Q_OBJECT

public:
    explicit ModelWidget2(QWidget *parent = nullptr);
    ~ModelWidget2();

    ModelCurveData calculateTheoreticalCurve(const QMap<QString, double>& params,
                                             const QVector<double>& providedTime = QVector<double>());

    void setHighPrecision(bool high);

private slots:
    void onCalculateClicked();
    void onResetParameters();
    void onExportData();
    void onExportImage();
    void onResetView();
    void onFitToData();
    void onChartSettings();
    void onDependentParamsChanged();
    void onShowPointsToggled(bool checked);

signals:
    void calculationCompleted(const QString &analysisType, const QMap<QString, double> &results);

private:
    void initChart();
    void setupConnections();

    QVector<double> parseInput(const QString& text);
    void setInputText(QLineEdit* edit, double value);
    void runCalculation();
    void plotCurve(const ModelCurveData& data, const QString& name, QColor color, bool isSensitivity);

    void calculatePDandDeriv(const QVector<double>& tD, const QMap<QString, double>& params,
                             std::function<double(double, const QMap<QString, double>&)> laplaceFunc,
                             QVector<double>& outPD, QVector<double>& outDeriv);

    double flaplace_composite(double z, const QMap<QString, double>& p);
    double PWD_inf(double z, double fs1, double fs2, double M12, double LfD, double rmD, int nf, const QVector<double>& xwD);
    double scaled_besseli(int v, double x);
    double adaptiveGauss(std::function<double(double)> f, double a, double b, double eps, int depth, int maxDepth);
    double gauss15(std::function<double(double)> f, double a, double b);
    double stefestCoefficient(int i, int N);
    double factorial(int n);

private:
    Ui::ModelWidget2 *ui;
    MouseZoom *m_plot;
    QCPTextElement *m_plotTitle;
    bool m_highPrecision;

    QVector<double> res_tD;
    QVector<double> res_pD;
    QVector<double> res_dpD;

    QList<QColor> m_colorList;
};

#endif // MODELWIDGET2_H
