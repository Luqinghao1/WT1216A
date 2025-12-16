#include "modelwidget4.h"
#include "ui_modelwidget4.h"

ModelWidget4::ModelWidget4(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::ModelWidget4)
{
    ui->setupUi(this);
}

ModelWidget4::~ModelWidget4()
{
    delete ui;
}

ModelCurveData ModelWidget4::calculateTheoreticalCurve(const QMap<QString, double>& params,
                                                       const QVector<double>& providedTime)
{
    return std::make_tuple(QVector<double>(), QVector<double>(), QVector<double>());
}
