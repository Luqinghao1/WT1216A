#ifndef MODELWIDGET5_H
#define MODELWIDGET5_H

#include <QWidget>
#include <QMap>
#include <QVector>
#include <tuple>

typedef std::tuple<QVector<double>, QVector<double>, QVector<double>> ModelCurveData;

namespace Ui {
class ModelWidget5;
}

class ModelWidget5 : public QWidget
{
    Q_OBJECT

public:
    explicit ModelWidget5(QWidget *parent = nullptr);
    ~ModelWidget5();

    ModelCurveData calculateTheoreticalCurve(const QMap<QString, double>& params,
                                             const QVector<double>& providedTime = QVector<double>());

private:
    Ui::ModelWidget5 *ui;
};

#endif // MODELWIDGET5_H
