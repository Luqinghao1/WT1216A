#ifndef MODELWIDGET6_H
#define MODELWIDGET6_H

#include <QWidget>
#include <QMap>
#include <QVector>
#include <tuple>

typedef std::tuple<QVector<double>, QVector<double>, QVector<double>> ModelCurveData;

namespace Ui {
class ModelWidget6;
}

class ModelWidget6 : public QWidget
{
    Q_OBJECT

public:
    explicit ModelWidget6(QWidget *parent = nullptr);
    ~ModelWidget6();

    ModelCurveData calculateTheoreticalCurve(const QMap<QString, double>& params,
                                             const QVector<double>& providedTime = QVector<double>());

private:
    Ui::ModelWidget6 *ui;
};

#endif // MODELWIDGET6_H
