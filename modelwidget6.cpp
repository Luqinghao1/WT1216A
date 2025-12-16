#include "modelwidget6.h"
#include "ui_modelwidget6.h"

ModelWidget6::ModelWidget6(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::ModelWidget6)
{
    ui->setupUi(this);
}

ModelWidget6::~ModelWidget6()
{
    delete ui;
}

ModelCurveData ModelWidget6::calculateTheoreticalCurve(const QMap<QString, double>& params,
                                                       const QVector<double>& providedTime)
{
    return std::make_tuple(QVector<double>(), QVector<double>(), QVector<double>());
}
