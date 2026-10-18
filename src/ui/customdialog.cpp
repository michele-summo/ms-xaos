#include <QtWidgets>
#define __USE_MINGW_ANSI_STDIO 1 // for long double support on Windows
#include <cstdio>

#include "customdialog.h"

#include "config.h"
#include "ui.h"
#include "misc-f.h"
#include "filter.h"
#include "i18n.h"
#include "ui_helper.h"
#include "viewtools.h"

#ifdef USE_FLOAT128
#include <quadmath.h>
#endif

QStringList fnames = {};

/* Where the exponent form starts to be the readable one.
 *
 * %g decides this for itself, by comparing the exponent against the number of
 * significant digits asked for -- which is no use when the digits are being
 * chosen for shortness: twenty came out as 2e+01, since one digit is all it
 * takes to read back as twenty. So the notation is chosen by size, and the
 * digits within it by what reads back. */
#define PLAIN_ABOVE ((number_t)1e12)  /* a thousand billion */
#define PLAIN_BELOW ((number_t)1e-7)  /* a hundred billionths */

static bool wants_exponent(number_t number)
{
    number_t size = number < 0 ? -number : number;
    if (size == 0)
        return false;
    return size >= PLAIN_ABOVE || size < PLAIN_BELOW;
}

/* One printing: `digits` significant ones in exponent form, or that many
 * places after the point in plain form. */
static void print_digits(char *buf, size_t size, number_t number, int digits,
                         bool exponent)
{
    char fs[16];
    const char form = exponent ? 'e' : 'f';
#ifdef USE_FLOAT128
    snprintf(fs, sizeof(fs), "%%.%iQ%c", digits, form);
    quadmath_snprintf(buf, size, fs, (__float128)number);
#else
#ifdef USE_LONG_DOUBLE
    snprintf(fs, sizeof(fs), "%%.%iL%c", digits, form);
    snprintf(buf, size, fs, (long double)number);
#else
    snprintf(fs, sizeof(fs), "%%.%i%c", digits, form);
    snprintf(buf, size, fs, (double)number);
#endif
#endif
}

QString CustomDialog::format(number_t number)
{
    char buf[256];
    char *end;

    /* The shortest text that reads back as this very number.
     *
     * Every digit the build carries used to be printed, so that a field could
     * not lose precision when it was read back. That is right for a
     * coordinate at a deep zoom and noise everywhere else: a convergence of
     * 1e-18 came back looking like 9.99999999999999999978e-19 and a
     * thousandth like 0.000999999999999999999958, because the decimal that
     * was typed is not one a binary float holds exactly and twenty-one digits
     * say so at length.
     *
     * Trying the lengths in turn and stopping at the first that reads back
     * unchanged keeps the guarantee exactly -- that is what reading back
     * unchanged means -- while saying 1e-18 where 1e-18 is what it holds, and
     * still spelling out every digit of a coordinate that needs them. */
    const bool exponent = wants_exponent(number);
    /* Plain form counts places after the point, so a small number needs more
     * of them than it has significant digits; exponent form counts the digits
     * themselves and never needs more than the build carries. */
    const int most = exponent ? NUMBER_DIGITS : NUMBER_DIGITS + 12;
    for (int digits = 0; digits < most; digits++) {
        print_digits(buf, sizeof(buf), number, digits, exponent);
        if (xstrtonum(buf, &end) == number && end != buf && *end == 0)
            return QString(buf);
    }
    print_digits(buf, sizeof(buf), number, most, exponent);
    return QString(buf);
}

/* The name the widget of field i goes by, which accept() finds it by again.
 * It was the field's label, so two fields with the same label -- the ITER of
 * two groups -- were one name, and both read back from the first. */
static QString fieldName(int i) { return QString("field%1").arg(i); }

CustomDialog::CustomDialog(struct uih_context *uih, const menuitem *item,
                           const menudialog *dialog, QWidget *parent)
    : QDialog(parent)
{
    m_menuitem = item;
    m_dialog = dialog;
    m_parameters = 0;
    turning = false;

    setWindowTitle(item->name);

    QBoxLayout *dialogLayout = new QBoxLayout(QBoxLayout::TopToBottom, this);
    /* The columns side by side, each a form held to the top; most dialogs
     * have the one. See menudialog's column. */
    QBoxLayout *columnsLayout = new QBoxLayout(QBoxLayout::LeftToRight);
    QFormLayout *formLayout = NULL;
    /* The form each field went into, for the fields shown only with a
     * choice, which are hidden and shown in their own form. */
    QList<QFormLayout *> fieldForm;
    auto addColumn = [this, columnsLayout, &formLayout]() {
        if (formLayout != NULL) {
            QFrame *rule = new QFrame(this);
            rule->setFrameShape(QFrame::VLine);
            rule->setFrameShadow(QFrame::Sunken);
            columnsLayout->addWidget(rule);
        }
        formLayout = new QFormLayout();
        QBoxLayout *column = new QBoxLayout(QBoxLayout::TopToBottom);
        column->addLayout(formLayout);
        column->addStretch();
        columnsLayout->addLayout(column);
    };
    addColumn();

    for (int i = 0; dialog[i].question; i++) {
        QString label(dialog[i].question);
        const QString name = fieldName(i);
        if (dialog[i].column && i > 0)
            addColumn();
        fieldForm.append(formLayout);
        if (dialog[i].separated) {
            if (!dialog[i].column) {
                QFrame *rule = new QFrame(this);
                rule->setFrameShape(QFrame::HLine);
                rule->setFrameShadow(QFrame::Sunken);
                formLayout->addRow(rule);
            }
            if (dialog[i].heading) {
                QLabel *heading = new QLabel(dialog[i].heading, this);
                QFont bold = heading->font();
                bold.setBold(true);
                heading->setFont(bold);
                formLayout->addRow(heading);
            }
        }
        if (dialog[i].type == DIALOG_COORD) {

            QLineEdit *real = new QLineEdit(format(dialog[i].deffloat), this);
            QFontMetrics metric(real->font());
            real->setMinimumWidth(metric.horizontalAdvance(real->text()) * 1.1);
            real->setObjectName(name + "real");
            // real->setValidator(new QDoubleValidator(real));

            QLineEdit *imag = new QLineEdit(format(dialog[i].deffloat2), this);
            imag->setObjectName(name + "imag");
            imag->setMinimumWidth(metric.horizontalAdvance(imag->text()) * 1.1);
            // imag->setValidator(new QDoubleValidator(imag));

            QBoxLayout *layout = new QBoxLayout(QBoxLayout::LeftToRight);
            layout->setContentsMargins(0, 0, 0, 0);
            layout->addWidget(real);
            layout->addWidget(new QLabel("+", this));
            layout->addWidget(imag);
            layout->addWidget(new QLabel("i", this));

            formLayout->addRow(label, layout);

        } else if (dialog[i].type == DIALOG_IFILE ||
                   dialog[i].type == DIALOG_OFILE) {

            QLineEdit *filename = new QLineEdit(dialog[i].defstr, this);
            QFontMetrics metric(filename->font());
            filename->setMinimumWidth(metric.horizontalAdvance(filename->text()) * 1.1);
            filename->setObjectName(name);

            QToolButton *chooser = new QToolButton(this);
            chooser->setProperty("field", i);
            chooser->setText("...");

            if (dialog[i].type == DIALOG_IFILE)
                connect(chooser, SIGNAL(clicked()), this,
                        SLOT(chooseInputFile()));
            else
                connect(chooser, SIGNAL(clicked()), this,
                        SLOT(chooseOutputFile()));

            QBoxLayout *layout = new QBoxLayout(QBoxLayout::LeftToRight);
            layout->setContentsMargins(0, 0, 0, 0);
            layout->addWidget(filename);
            layout->addWidget(chooser);

            formLayout->addRow(label, layout);

        } else if (dialog[i].type == DIALOG_IFILES) {

            QTextEdit *filenames = new QTextEdit(this);
            filenames->setMaximumHeight(filenames->height() * 2);
            filenames->setObjectName(name);

            QToolButton *chooser = new QToolButton(this);
            chooser->setProperty("field", i);
            chooser->setText("...");
            connect(chooser, SIGNAL(clicked()), this,
                    SLOT(chooseInputFiles()));

            QBoxLayout *layout = new QBoxLayout(QBoxLayout::LeftToRight);
            layout->setContentsMargins(0, 0, 0, 0);
            layout->addWidget(filenames);
            layout->addWidget(chooser);

            formLayout->addRow(label, layout);

        } else if (dialog[i].type == DIALOG_CHOICE) {

            QComboBox *combo = new QComboBox(this);
            combo->setObjectName(name);

            const char **str = (const char **)dialog[i].defstr;
            for (int j = 0; str[j] != NULL; j++)
                combo->addItem(TR("Dialog", str[j]));
            combo->setCurrentIndex(dialog[i].defint);

            formLayout->addRow(label, combo);

        } else if (dialog[i].type == DIALOG_PALSLIDER) {

            gradientpal = clonepalette(uih->image->palette);
            palcontext = uih;
            // 3 inputs decide color, Algorithm Number, Seed and shift
            // For Algorithm number
            algono = new QSpinBox(this);
            algono->setObjectName(label + "algono");
            algono->setValue(palcontext->palettetype);
            algono->setRange(1, PALGORITHMS);
            /* The number alone says nothing once there are more than the three
             * XaoS always had, so the name goes after it -- as the box's own
             * suffix, since the preview below takes its width from the box. */
            algono->setSuffix(QString("   ") +
                              TR("Palette", palette_algorithm_name(algono->value())));
            connect(algono, &QSpinBox::valueChanged, algono, [this](int n) {
                algono->setSuffix(QString("   ") +
                                  TR("Palette", palette_algorithm_name(n)));
            });

            // Algo Slider
            algoslider = new QSlider(Qt::Horizontal, this);
            algoslider->setObjectName(label);
            algoslider->setRange(1, PALGORITHMS);
            algoslider->setValue(algono->value());
            // algoslider->setMinimumWidth(this->width()*2);

            // For Seed Number
            seedno = new QSpinBox(this);
            seedno->setObjectName(label + "seedno");
            seedno->setRange(0, gradientpal->size);
            seedno->setValue(palcontext->paletteseed);

            // Seed Slider
            seedslider = new QSlider(Qt::Horizontal, this);
            seedslider->setObjectName(label);
            seedslider->setRange(0, gradientpal->size);
            seedslider->setValue(seedno->value());

            // For Shift Number
            shiftno = new QSpinBox(this);
            shiftno->setObjectName(label + "shiftno");
            shiftno->setRange(0, gradientpal->size);
            shiftno->setValue(palcontext->paletteshift + palcontext->manualpaletteshift);

            // Shift Slider
            shiftslider = new QSlider(Qt::Horizontal, this);
            shiftslider->setObjectName(label);
            shiftslider->setRange(0, gradientpal->size);
            shiftslider->setValue(shiftno->value());

            // Add them to Layout
            formLayout->addRow(TR("Dialog", "Algorithm"), algono);
            formLayout->addWidget(algoslider);
            formLayout->addRow(TR("Dialog", "Seed"), seedno);
            formLayout->addWidget(seedslider);
            formLayout->addRow(TR("Dialog", "Shift"), shiftno);
            formLayout->addWidget(shiftslider);

            img = new QLabel(this);
            img->setScaledContents(true);
            formLayout->addRow(img);
            updateVisualiser();
            connect(algono,SIGNAL(valueChanged(int)), algoslider, SLOT(setValue(int)));
            connect(algoslider, SIGNAL(valueChanged(int)), algono, SLOT(setValue(int)));
            connect(algono, SIGNAL(valueChanged(int)), this, SLOT(updateVisualiser()));
            connect(seedno,SIGNAL(valueChanged(int)), seedslider, SLOT(setValue(int)));
            connect(seedslider, SIGNAL(valueChanged(int)), seedno, SLOT(setValue(int)));
            connect(seedno, SIGNAL(valueChanged(int)), this, SLOT(updateVisualiser()));
            connect(shiftno,SIGNAL(valueChanged(int)), shiftslider, SLOT(setValue(int)));
            connect(shiftslider, SIGNAL(valueChanged(int)), shiftno, SLOT(setValue(int)));
            connect(shiftno, SIGNAL(valueChanged(int)), this, SLOT(updateVisualiser()));

        } else if (dialog[i].type == DIALOG_PALPICKER) {

            palcontext = uih;
            getDEFSEGMENTColor(newColors);

            QList< QPushButton* > buttons;
            QBoxLayout *layout1 = new QBoxLayout(QBoxLayout::LeftToRight);
            QBoxLayout *layout2 = new QBoxLayout(QBoxLayout::LeftToRight);
            QBoxLayout *layout3 = new QBoxLayout(QBoxLayout::LeftToRight);
            for(auto bidx = 0; bidx < 31; ++bidx ) {
                auto button = new QPushButton{ QString::number(bidx) };
                button->setObjectName(QString::number(bidx));
                QColor color(newColors[bidx][0], newColors[bidx][1], newColors[bidx][2]);
                QPalette pal = button->palette();
                button->setAutoFillBackground(true);
                pal.setColor(QPalette::Button, color);
                button->setPalette(pal);
                button->update();
                buttons << button;
                if(bidx <= 10)
                    layout1->addWidget(button);
                else if(bidx>10 and bidx <= 20)
                    layout2->addWidget(button);
                else
                    layout3->addWidget(button);

                connect(button, SIGNAL(clicked()), this, SLOT(colorPicker()));
            }
            formLayout->addRow(layout1);
            formLayout->addRow(layout2);
            formLayout->addRow(layout3);

        } else if (dialog[i].type == DIALOG_LIST) {

            QComboBox *list = new QComboBox(this);
            list->setObjectName(name);
            list->setEditable(true);
            list->addItem(dialog[i].defstr);

            QSettings settings;
            QStringList formulas = settings.value("Formulas/UserFormulas").toStringList();
            list->addItems(formulas);

            formLayout->addRow(label, list);

        } else {

            QLineEdit *field = new QLineEdit(this);
            field->setObjectName(name);

            if (dialog[i].type == DIALOG_INT) {
                field->setText(QString::number(dialog[i].defint));
                field->setValidator(new QIntValidator(field));
            } else if (dialog[i].type == DIALOG_FLOAT) {
                field->setText(format(dialog[i].deffloat));
                // field->setValidator(new QDoubleValidator(field));
            } else {
                field->setText(dialog[i].defstr);
            }
            QFontMetrics metric(field->font());
            field->setMinimumWidth(metric.horizontalAdvance(field->text()) * 1.1);
            formLayout->addRow(label, field);
        }

        if (dialog[i].viewtools)
            addViewTools(formLayout, i, i - dialog[i].viewtools);
    }

    /* The fields shown only while a choice before them holds a given entry --
     * the radius of a group, in smooth alone -- following the choice as it
     * is changed. A hidden field is still read back on OK, with what it
     * held. */
    for (int i = 0; dialog[i].question; i++) {
        const int back = dialog[i].shownwith;
        if (back <= 0 || back > i)
            continue;
        QComboBox *choice = findChild<QComboBox *>(fieldName(i - back));
        QWidget *field = findChild<QWidget *>(fieldName(i));
        if (choice == NULL || field == NULL)
            continue;
        const int at = dialog[i].shownat;
        QFormLayout *form = fieldForm[i];
        auto follow = [this, form, field, at](int index) {
            form->setRowVisible(field, index == at);
            adjustSize();
        };
        follow(choice->currentIndex());
        connect(choice, &QComboBox::currentIndexChanged, this, follow);
    }

    dialogLayout->addLayout(columnsLayout);

    QDialogButtonBox *buttonBox =
        new QDialogButtonBox((QDialogButtonBox::Ok | QDialogButtonBox::Cancel),
                             Qt::Horizontal, this);
    connect(buttonBox, SIGNAL(accepted()), this, SLOT(accept()));
    connect(buttonBox, SIGNAL(rejected()), this, SLOT(reject()));
    dialogLayout->addWidget(buttonBox);

    setLayout(dialogLayout);
}

void CustomDialog::accept()
{
    int nitems;
    char *ps;
    for (nitems = 0; m_dialog[nitems].question; nitems++)
        ;
    m_parameters = (dialogparam *)malloc(sizeof(*m_parameters) * nitems);

    for (int i = 0; i < nitems; i++) {
        const QString name = fieldName(i);

        if (m_dialog[i].type == DIALOG_COORD) {

            QLineEdit *real = findChild<QLineEdit *>(name + "real");
            QLineEdit *imag = findChild<QLineEdit *>(name + "imag");

            m_parameters[i].dcoord[0] = xstrtonum(real->text().toUtf8(), &ps);
            m_parameters[i].dcoord[1] = xstrtonum(imag->text().toUtf8(), &ps);

        } else if (m_dialog[i].type == DIALOG_CHOICE) {

            QComboBox *field = findChild<QComboBox *>(name);
            m_parameters[i].dint = field->currentIndex();

        }
        else if (m_dialog[i].type == DIALOG_IFILES){

            QTextEdit *field = findChild<QTextEdit *>(name);
            QString raw_fnames = field->toPlainText();
            fnames = raw_fnames.split("\n");

        } else {

            QLineEdit *field = findChild<QLineEdit *>(name);

            if (m_dialog[i].type == DIALOG_INT)
                m_parameters[i].dint = field->text().toInt();
            else if (m_dialog[i].type == DIALOG_FLOAT)
                m_parameters[i].number = xstrtonum(field->text().toUtf8(), &ps);
            else if (m_dialog[i].type == DIALOG_PALSLIDER) {
                // Copy data to original context
                palcontext->manualpaletteshift = 0;
                palcontext->palettetype = algono->value();
                palcontext->paletteseed = seedno->value();
                palcontext->paletteshift = shiftno->value();
                m_parameters[i].dint = 1;
                destroypalette(gradientpal);
            } else if (m_dialog[i].type == DIALOG_PALPICKER) {
                mkcustompalette(palcontext->image->palette, newColors);
            } else if (m_dialog[i].type == DIALOG_LIST) {
                QComboBox *list = findChild<QComboBox *>(name);
                m_parameters[i].dstring = strdup(list->currentText().toUtf8());
            }
            else
                m_parameters[i].dstring = strdup(field->text().toUtf8());
        }
    }

    QDialog::accept();
}

dialogparam *CustomDialog::parameters() { return m_parameters; }

/* The slices a full turn is cut into, kept from one opening of the View
 * dialog to the next while XaoS runs. Nowhere else: it is not part of the
 * view, and a saved position says nothing of it. */
static int viewSlices = 0;

/* What a field holds, or false where it does not read as a number -- so that
 * a tool leaves alone what it cannot make sense of. */
static bool readNumber(QLineEdit *field, number_t *number)
{
    const QByteArray text = field->text().trimmed().toUtf8();
    char *end;
    *number = xstrtonum(text.constData(), &end);
    return end != text.constData() && *end == 0;
}

/* A rule and the slices under it, a rule and the axes under that: see
 * menudialog's viewtools. They change what the angle and center fields hold
 * and nothing else, so OK applies it all as one step and Cancel none of it. */
void CustomDialog::addViewTools(QFormLayout *form, int angle, int center)
{
    QLineEdit *angleField = findChild<QLineEdit *>(fieldName(angle));
    QLineEdit *real = findChild<QLineEdit *>(fieldName(center) + "real");
    QLineEdit *imag = findChild<QLineEdit *>(fieldName(center) + "imag");
    if (angleField == NULL || real == NULL || imag == NULL)
        return;

    QFrame *rule = new QFrame(this);
    rule->setFrameShape(QFrame::HLine);
    rule->setFrameShadow(QFrame::Sunken);
    form->addRow(rule);

    QSpinBox *slices = new QSpinBox(this);
    slices->setRange(0, 999);
    slices->setValue(viewSlices);
    slices->setToolTip(TR("Dialog", "How many equal slices a full turn is cut "
                                    "into: - and + turn the angle by one"));
    QToolButton *back = new QToolButton(this);
    back->setText(QString(QChar(0x2212)));
    back->setToolTip(TR("Dialog", "Turn the angle back by one slice, "
                                  "360 degrees over the slices"));
    QToolButton *on = new QToolButton(this);
    on->setText("+");
    on->setToolTip(TR("Dialog", "Turn the angle on by one slice, "
                                "360 degrees over the slices"));

    /* Nought slices is no slice to turn by. */
    auto enable = [back, on](int n) {
        back->setEnabled(n > 0);
        on->setEnabled(n > 0);
    };
    enable(slices->value());
    connect(slices, &QSpinBox::valueChanged, this, [this, enable](int n) {
        viewSlices = n;
        turning = false;
        enable(n);
    });
    /* Typed by hand, the angle is where turning starts from again. setText
     * below is not typing, so the slices turned stay counted from where
     * they began. */
    connect(angleField, &QLineEdit::textEdited, this,
            [this](const QString &) { turning = false; });
    auto turn = [this, angleField, slices](int step) {
        if (!turning) {
            if (!readNumber(angleField, &turnfrom))
                return;
            turnsteps = 0;
            turning = true;
        }
        turnsteps += step;
        angleField->setText(
            format(view_turn(turnfrom, turnsteps, slices->value())));
    };
    connect(back, &QToolButton::clicked, this, [turn] { turn(-1); });
    connect(on, &QToolButton::clicked, this, [turn] { turn(1); });

    QBoxLayout *turnRow = new QBoxLayout(QBoxLayout::LeftToRight);
    turnRow->setContentsMargins(0, 0, 0, 0);
    turnRow->addWidget(slices);
    turnRow->addWidget(back);
    turnRow->addWidget(on);
    turnRow->addStretch();
    form->addRow(TR("Dialog", "Slices:"), turnRow);

    rule = new QFrame(this);
    rule->setFrameShape(QFrame::HLine);
    rule->setFrameShadow(QFrame::Sunken);
    form->addRow(rule);

    QToolButton *alongR = new QToolButton(this);
    alongR->setText("R");
    alongR->setToolTip(TR("Dialog", "Move the view across the picture until "
                                    "the origin is in the middle from left to "
                                    "right, at the angle above"));
    QToolButton *alongI = new QToolButton(this);
    alongI->setText("I");
    alongI->setToolTip(TR("Dialog", "Move the view up or down the picture "
                                    "until the origin is in the middle from "
                                    "top to bottom, at the angle above"));
    auto centre = [angleField, real, imag](int axis) {
        number_t a, cr, ci;
        if (!readNumber(angleField, &a) || !readNumber(real, &cr) ||
            !readNumber(imag, &ci))
            return;
        view_centeronaxis(a, axis, &cr, &ci);
        real->setText(format(cr));
        imag->setText(format(ci));
    };
    connect(alongR, &QToolButton::clicked, this,
            [centre] { centre(VIEW_AXIS_R); });
    connect(alongI, &QToolButton::clicked, this,
            [centre] { centre(VIEW_AXIS_I); });

    QBoxLayout *axisRow = new QBoxLayout(QBoxLayout::LeftToRight);
    axisRow->setContentsMargins(0, 0, 0, 0);
    axisRow->addWidget(alongR);
    axisRow->addWidget(alongI);
    axisRow->addStretch();
    form->addRow(TR("Dialog", "Center by axis:"), axisRow);
}

/* The field a "..." button belongs to, and the question it asks. */
static int chooserField(QObject *chooser)
{
    return chooser->property("field").toInt();
}

void CustomDialog::chooseInputFile()
{
    const int i = chooserField(sender());
    QLineEdit *field = findChild<QLineEdit *>(fieldName(i));

    QSettings settings;
    QString fileLocation = settings.value("MainWindow/lastFileLocation", QDir::homePath()).toString();
    QString fileName = QFileDialog::getOpenFileName(
            this, m_dialog[i].question, fileLocation, "*.xpf *.png *.xaf");
    if (!fileName.isNull()) {
        field->setText(fileName);
        settings.setValue("MainWindow/lastFileLocation", QFileInfo(fileName).absolutePath());
    }
}

void CustomDialog::chooseOutputFile()
{
    const int i = chooserField(sender());
    QLineEdit *field = findChild<QLineEdit *>(fieldName(i));
    QSettings settings;
    QString fileLocation = settings.value("MainWindow/lastFileLocation", QDir::homePath()).toString();
    QString fileName = QFileDialog::getSaveFileName(
        this, m_dialog[i].question, fileLocation);
    if (!fileName.isNull()) {
        field->setText(fileName);
        settings.setValue("MainWindow/lastFileLocation", QFileInfo(fileName).absolutePath());
    }
}

void CustomDialog::chooseInputFiles()
{
    const int i = chooserField(sender());
    QTextEdit *field = findChild<QTextEdit *>(fieldName(i));
    QSettings settings;
    QString fileLocation = settings.value("MainWindow/lastFileLocation", QDir::homePath()).toString();
    QStringList fileNames = QFileDialog::getOpenFileNames(
        this, m_dialog[i].question, fileLocation, "*.xpf *.xaf");
    if(!fileNames.isEmpty()) {
        for(auto file: fileNames) {
            field->append(file);
        }
        settings.setValue("MainWindow/lastFileLocation", QFileInfo(fileNames[0]).absolutePath());
    }
}

void CustomDialog::updateVisualiser()
{
    // Get updated Colors
    int colors[101][3];
    getPaletteColor(gradientpal, seedno->value(),
                    algono->value()-1 < 0? 0:algono->value()-1, shiftno->value(), colors);

    // Load Curve
    QImage palImage(100, 1, QImage::Format_RGB32);

    // Fill Curve
    for(int i=0;i<100;i++) {
        QRgb value = qRgb(colors[i][0], colors[i][1], colors[i][2]);
        palImage.setPixelColor(i, 0, value);
    }

    // Save Result
    QPixmap newImage = QPixmap::fromImage(palImage.scaled(this->algono->width(),
                                                          this->algono->height()));
    img->setPixmap(newImage);
}

void CustomDialog::colorPicker()
{
    QPushButton* button = qobject_cast<QPushButton*>(sender());
    int idx = button->objectName().toInt();
    QColor color = QColorDialog::getColor(QColor(newColors[idx][0], newColors[idx][1],
            newColors[idx][2]), this);
    QPalette pal = button->palette();
    button->setAutoFillBackground(true);
    pal.setColor(QPalette::Button, color);
    button->setPalette(pal);
    button->update();
    newColors[idx][0] = color.red();
    newColors[idx][1] = color.green();
    newColors[idx][2] = color.blue();
}
