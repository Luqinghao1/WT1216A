#ifndef MODELWIDGET4_H
#define MODELWIDGET4_H

#include <QWidget>
#include <QMap>
#include <QVector>
#include <tuple>

typedef std::tuple<QVector<double>, QVector<double>, QVector<double>> ModelCurveData;

namespace Ui {
class ModelWidget4;
}

class ModelWidget4 : public QWidget
{
    Q_OBJECT

public:
    explicit ModelWidget4(QWidget *parent = nullptr);
    ~ModelWidget4();

    ModelCurveData calculateTheoreticalCurve(const QMap<QString, double>& params,
                                             const QVector<double>& providedTime = QVector<double>());

private:
    Ui::ModelWidget4 *ui;
};

#endif // MODELWIDGET4_H
