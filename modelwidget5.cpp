#include "modelwidget5.h"
#include "ui_modelwidget5.h"

ModelWidget5::ModelWidget5(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::ModelWidget5)
{
    ui->setupUi(this);
}

ModelWidget5::~ModelWidget5()
{
    delete ui;
}

ModelCurveData ModelWidget5::calculateTheoreticalCurve(const QMap<QString, double>& params,
                                                       const QVector<double>& providedTime)
{
    return std::make_tuple(QVector<double>(), QVector<double>(), QVector<double>());
}
