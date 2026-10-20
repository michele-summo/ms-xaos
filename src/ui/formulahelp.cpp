/* Help -> User formula reference: a window listing what a user-defined
 * formula may contain.
 *
 * A window of its own and not a modal dialog, so that it can be left open at
 * the side while a formula is being written -- which is the only way a
 * reference is any use. It is built once and re-shown afterwards, so reopening
 * it returns to the tab and the scroll position it was left at.
 *
 * The content comes from formulahelp_data.cpp, which is generated from the
 * parser's own function table and holds no Qt; see formulahelp.h.
 */

#include <QtWidgets>

#include "config.h"
#include "formulahelp.h"

/* The line under the formula bar: what the call the cursor is standing in
 * takes, with the argument being written picked out.
 *
 * The bar is one line and the reference is behind a menu, so a formula with
 * six optional arguments has to be written from memory. This says what the
 * places are as one moves between them. Which call and which argument is
 * formula_hint_at, which is written without Qt and tested; what is left here
 * is the wording, which is why it lives beside the reference it words rather
 * than beside the bar it is shown under.
 *
 * Rich text, since the emphasis is the whole point -- a signature with six
 * arguments says nothing about which of them is under the cursor. Everything
 * put in is escaped: the reference holds no markup today, and an argument that
 * grew a "<" would otherwise vanish instead of being shown.
 */
QString formula_hint_text(const QString &text, int cursor)
{
    /* The cursor counts in the units QString counts in, and the walk counts
     * bytes, so the prefix is measured rather than assumed equal. */
    const QByteArray utf8 = text.toUtf8();
    const int at = text.left(cursor).toUtf8().size();

    struct formula_hint hint = formula_hint_at(utf8.constData(), at);
    if (!hint.fn)
        return QString();

    const QString name = QString::fromUtf8(hint.fn->name).toHtmlEscaped();
    const char *args = hint.fn->args ? hint.fn->args : "";
    int start = 0, len = 0;
    if (hint.arg < 0 || !formula_hint_arg_span(args, hint.arg, &start, &len))
        start = len = 0;

    const QString all = QString::fromUtf8(args);
    const QString before = QString::fromUtf8(args, start).toHtmlEscaped();
    const QString here = QString::fromUtf8(args + start, len).toHtmlEscaped();
    const QString after = QString::fromUtf8(args + start + len).toHtmlEscaped();

    QString shown = len ? before + "<b>" + here + "</b>" + after
                        : all.toHtmlEscaped();
    return name + "(" + shown + ")";
}

/* A table that measures its rows again whenever its width changes.
 *
 * A wrapped cell is as tall as the width lets it be, and the width is not
 * known while the table is being built: resizeRowsToContents there measures
 * against a width the table will not have, so the first layout was right only
 * by accident and every resize afterwards left the descriptions either clipped
 * or floating in rows far too tall. Measuring on the resize is what makes the
 * columns fit the window.
 */
class WrappingTable : public QTableWidget
{
  public:
    WrappingTable(int rows, int columns) : QTableWidget(rows, columns) {}

  protected:
    void resizeEvent(QResizeEvent *event) override
    {
        QTableWidget::resizeEvent(event);
        resizeRowsToContents();
    }
};

/* One tab: a two-column table of names and descriptions, with the section
 * rows spanning both columns.
 *
 * QTableWidget rather than anything lighter because it scrolls by itself and
 * keeps the header in place while it does, which is what makes a list of
 * ninety-six entries usable. */
static QWidget *buildTable(const struct formula_help_row *rows,
                           const QString &nameHeading, bool withArgs = false)
{
    int count = 0;
    for (const struct formula_help_row *r = rows; r->name || r->section; r++)
        count++;

    const int columns = withArgs ? 3 : 2;
    QTableWidget *table = new WrappingTable(count, columns);
    QStringList headings;
    headings << nameHeading;
    if (withArgs)
        headings << QObject::tr("Takes");
    headings << QObject::tr("What it does");
    table->setHorizontalHeaderLabels(headings);
    table->verticalHeader()->setVisible(false);
    table->setEditTriggers(QAbstractItemView::NoEditTriggers);
    table->setSelectionBehavior(QAbstractItemView::SelectRows);
    table->setAlternatingRowColors(true);
    table->setWordWrap(true);
    /* Nothing is cut short and nothing is pushed sideways: the description
     * column wraps to whatever width is left instead. */
    table->setTextElideMode(Qt::ElideNone);
    table->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    /* The name column is as wide as the longest name and no wider.
     *
     * Sizing it to its contents instead is what burst the window: a section
     * heading is a cell of column zero spanning into column one, so the
     * heading counted as content of the name column and made it as wide as
     * the longest sentence in the table, leaving the descriptions no room at
     * all -- the second column simply disappeared. Measuring the names is
     * both what was meant and immune to whatever the headings say. */
    QFont fixed = QFontDatabase::systemFont(QFontDatabase::FixedFont);
    QFontMetrics metrics(fixed);
    int names = metrics.horizontalAdvance(nameHeading);
    for (const struct formula_help_row *r = rows; r->name || r->section; r++)
        if (r->name)
            names = qMax(names, metrics.horizontalAdvance(QString(r->name)));
    table->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Fixed);
    table->setColumnWidth(0, names + 3 * metrics.horizontalAdvance(QChar('m')));
    if (withArgs) {
        /* Measured the same way and for the same reason as the names. */
        int taken = metrics.horizontalAdvance(QObject::tr("Takes"));
        for (const struct formula_help_row *r = rows; r->name || r->section; r++)
            if (r->args)
                taken = qMax(taken, metrics.horizontalAdvance(QString(r->args)));
        /* Capped, and wrapping past the cap: a call with five optional
         * arguments and their defaults written out is longer than any window
         * should give to a column that is not the description. */
        int cap = 34 * metrics.horizontalAdvance(QChar('m'));
        table->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Fixed);
        table->setColumnWidth(
            1, qMin(taken + 3 * metrics.horizontalAdvance(QChar('m')), cap));
    }
    table->horizontalHeader()->setSectionResizeMode(columns - 1,
                                                    QHeaderView::Stretch);

    int row = 0;
    for (const struct formula_help_row *r = rows; r->name || r->section;
         r++, row++) {
        if (r->section != NULL) {
            QTableWidgetItem *item =
                new QTableWidgetItem(QCoreApplication::translate("FormulaHelp", r->section));
            QFont bold = item->font();
            bold.setBold(true);
            item->setFont(bold);
            item->setBackground(table->palette().alternateBase());
            table->setItem(row, 0, item);
            /* The heading is one cell across the lot. */
            table->setSpan(row, 0, 1, columns);
            continue;
        }
        QTableWidgetItem *name = new QTableWidgetItem(QString(r->name));
        name->setFont(fixed);
        table->setItem(row, 0, name);
        if (withArgs) {
            QTableWidgetItem *args = new QTableWidgetItem(QString(r->args));
            args->setFont(fixed);
            table->setItem(row, 1, args);
        }
        table->setItem(row, columns - 1,
                       new QTableWidgetItem(QCoreApplication::translate("FormulaHelp", r->summary)));
    }
    table->resizeRowsToContents();
    return table;
}

/* The Tilings tab: a picture of every tiling randsctile takes, captioned as
 * the Values tab lists it.
 *
 * Forty-five numbers with a line of words each say little about what the
 * shapes are, and a picture says it at a glance. The pictures are drawn by
 * randsctile itself, once, by tools/randsctile-thumbnails.cpp, and kept in
 * the binary, since they never change between one opening and the next;
 * formula_help_test checks they are still what it draws.
 *
 * A double click copies the start of the call, so that the tiling picked
 * goes into the formula without its number being looked up again. */
static QWidget *buildTilings()
{
    const struct formula_help_row *rows;
    const int count = formula_help_tilings(&rows);
    const int side = 128; /* the pictures are stored at twice this */

    QLabel *hint = new QLabel(QObject::tr(
        "Double-click a tiling to copy the start of its call to the "
        "clipboard."));
    hint->setTextFormat(Qt::PlainText);
    hint->setWordWrap(true);
    hint->setMargin(4);

    QListWidget *list = new QListWidget;
    list->setViewMode(QListView::IconMode);
    list->setMovement(QListView::Static);
    list->setResizeMode(QListView::Adjust);
    list->setWrapping(true);
    /* Not uniform sizes: those measure the first item for all of them, and
     * "1  squares" is one line, so every longer caption was cut to one. */
    list->setWordWrap(true);
    list->setSelectionMode(QAbstractItemView::SingleSelection);
    list->setIconSize(QSize(side, side));
    /* A cell as wide as a picture and a margin, and as tall as one and four
     * lines of caption: the longest takes three at this width. */
    QFontMetrics metrics(list->font());
    list->setGridSize(
        QSize(side + 3 * metrics.horizontalAdvance(QChar('m')),
              side + 4 * metrics.lineSpacing() + metrics.height() / 2));

    for (int k = 1; k <= count; k++) {
        QPixmap picture(
            QString(":/images/tilings/tiling-%1.png").arg(k, 2, 10, QChar('0')));
        picture.setDevicePixelRatio(2);
        /* The same picture when selected: left to itself Qt tints a selected
         * icon with the highlight, and the colours of the tiling with it. */
        QIcon icon(picture);
        icon.addPixmap(picture, QIcon::Selected);
        const QString what =
            QCoreApplication::translate("FormulaHelp", rows[k - 1].summary);
        QListWidgetItem *item = new QListWidgetItem(
            icon, QString("%1  %2").arg(k).arg(what), list);
        item->setData(Qt::UserRole, k);
        item->setToolTip(QString("randsctile(%1, seed, ...)\n%2").arg(k).arg(what));
    }

    QObject::connect(
        list, &QListWidget::itemDoubleClicked, hint,
        [hint](QListWidgetItem *item) {
            const QString call =
                QString("randsctile(%1, ").arg(item->data(Qt::UserRole).toInt());
            QGuiApplication::clipboard()->setText(call);
            hint->setText(QObject::tr("Copied \"%1\" to the clipboard: paste "
                                      "it into the formula and write the seed "
                                      "after it.")
                              .arg(call));
        });

    QWidget *page = new QWidget;
    QVBoxLayout *layout = new QVBoxLayout(page);
    layout->setContentsMargins(0, 0, 0, 0);
    layout->addWidget(hint);
    layout->addWidget(list);
    return page;
}

/* Why a field gives no number, worded for the line under the fields. Empty
 * is not among them: an empty length is where the tab starts, answered by
 * saying what to write, and an empty magnitude stands for its default. */
static QString fieldProblem(int reading)
{
    switch (reading) {
        case FORMULA_POLY_NOT_NUMBER:
            return QObject::tr("not a number");
        case FORMULA_POLY_NOT_WHOLE:
            return QObject::tr("not a whole number");
        case FORMULA_POLY_TOO_SHORT:
            return QObject::tr("at least 1");
        case FORMULA_POLY_TOO_LONG:
            return QObject::tr("at most %1: poly takes %2 arguments, and z is "
                               "one of them")
                .arg(FORMULA_POLY_LONGEST)
                .arg(FORMULA_POLY_LONGEST + 1);
        case FORMULA_POLY_NEGATIVE:
            return QObject::tr("a magnitude is not below 0");
        case FORMULA_POLY_TOO_LARGE:
            return QObject::tr("at most %1").arg(FORMULA_POLY_LARGEST);
    }
    return QString();
}

/* The Utilities tab: things that write a formula rather than describe one.
 *
 * For now a random poly: two lengths, two magnitudes and a delta of the
 * arguments, and a formula drawn from them whenever all five read -- an
 * empty delta reading as none -- drawn again on every change to any of
 * them, and on the button. While a length is empty or anything is not a
 * number there is no formula at all, not the last one left standing, so that
 * what is shown always answers what is written. The drawing is
 * formula_poly_random, without Qt and tested against the parser; this is only
 * the fields around it. */
static QWidget *buildUtilities()
{
    QFont fixed = QFontDatabase::systemFont(QFontDatabase::FixedFont);

    QGroupBox *box = new QGroupBox(QObject::tr("Random poly formula"));
    QLabel *about = new QLabel(QObject::tr(
        "A random polynomial in z, written as poly(z, k1, ..., km). Its length "
        "is m, the number of coefficients, and its degree m - 1. Each "
        "coefficient is a complex number with a magnitude drawn between the "
        "minimum and the maximum -- 0 and 1 when left empty -- and an argument "
        "drawn at random: round the whole turn, or, given a delta in degrees, "
        "the first round the whole turn and every other within delta degrees "
        "across, centred on the first's. There is no c in it: z starts at the pixel, so the "
        "picture is the polynomial's own filled Julia set, the same in "
        "mandelbrot and julia mode. Add c where you want it, as in "
        "poly(...) + c."));
    about->setTextFormat(Qt::PlainText);
    about->setWordWrap(true);

    QLineEdit *shortest = new QLineEdit;
    QLineEdit *longest = new QLineEdit;
    QLineEdit *smallest = new QLineEdit;
    QLineEdit *largest = new QLineEdit;
    QLineEdit *delta = new QLineEdit;
    shortest->setPlaceholderText(QObject::tr("e.g. %1").arg(3));
    longest->setPlaceholderText(QObject::tr("e.g. %1").arg(6));
    /* What an empty magnitude stands for, shown where it will be read. */
    smallest->setPlaceholderText("0");
    largest->setPlaceholderText("1");
    delta->setPlaceholderText(QObject::tr("empty: any argument"));
    QFormLayout *fields = new QFormLayout;
    fields->addRow(QObject::tr("Minimum length:"), shortest);
    fields->addRow(QObject::tr("Maximum length:"), longest);
    fields->addRow(QObject::tr("Minimum magnitude:"), smallest);
    fields->addRow(QObject::tr("Maximum magnitude:"), largest);
    fields->addRow(QObject::tr("Argument delta (degrees):"), delta);

    QLineEdit *formula = new QLineEdit;
    formula->setReadOnly(true);
    formula->setFont(fixed);
    QLabel *status = new QLabel;
    status->setTextFormat(Qt::PlainText);
    status->setWordWrap(true);

    QPushButton *again = new QPushButton(QObject::tr("Generate another"));
    QPushButton *copy = new QPushButton(QObject::tr("Copy"));
    QHBoxLayout *buttons = new QHBoxLayout;
    buttons->addWidget(again);
    buttons->addWidget(copy);
    buttons->addStretch();

    QVBoxLayout *inside = new QVBoxLayout(box);
    inside->addWidget(about);
    inside->addLayout(fields);
    inside->addWidget(formula);
    inside->addWidget(status);
    inside->addLayout(buttons);

    auto generate = [=]() {
        int lo = 0, hi = 0;
        double low = 0, high = 1, across = -1;
        const int a =
            formula_poly_length(shortest->text().toUtf8().constData(), &lo);
        const int b =
            formula_poly_length(longest->text().toUtf8().constData(), &hi);
        const int c = formula_poly_magnitude(
            smallest->text().toUtf8().constData(), 0, &low);
        const int d = formula_poly_magnitude(
            largest->text().toUtf8().constData(), 1, &high);
        const int e =
            formula_poly_delta(delta->text().toUtf8().constData(), &across);
        QString why;
        if (a == FORMULA_POLY_EMPTY || b == FORMULA_POLY_EMPTY)
            why = QObject::tr("Write the minimum and the maximum length, as "
                              "whole numbers from 1 to %1.")
                      .arg(FORMULA_POLY_LONGEST);
        else if (a != FORMULA_POLY_OK)
            why = QObject::tr("Minimum length: %1.").arg(fieldProblem(a));
        else if (b != FORMULA_POLY_OK)
            why = QObject::tr("Maximum length: %1.").arg(fieldProblem(b));
        else if (lo > hi)
            why = QObject::tr("The minimum length is more than the maximum.");
        else if (c != FORMULA_POLY_OK)
            why = QObject::tr("Minimum magnitude: %1.").arg(fieldProblem(c));
        else if (d != FORMULA_POLY_OK)
            why = QObject::tr("Maximum magnitude: %1.").arg(fieldProblem(d));
        else if (low > high)
            why = largest->text().trimmed().isEmpty()
                      ? QObject::tr("The minimum magnitude is more than the "
                                    "maximum, which is 1 when left empty.")
                      : QObject::tr("The minimum magnitude is more than the "
                                    "maximum.");
        else if (high == 0)
            why = QObject::tr("The maximum magnitude is 0: every coefficient "
                              "would be nought.");
        else if (e == FORMULA_POLY_NEGATIVE)
            why = QObject::tr("Argument delta: not below 0.");
        else if (e == FORMULA_POLY_TOO_WIDE)
            why = QObject::tr("Argument delta: at most %1, a whole turn.")
                      .arg(FORMULA_POLY_WIDEST);
        else if (e != FORMULA_POLY_OK)
            why = QObject::tr("Argument delta: %1.").arg(fieldProblem(e));
        const bool made = why.isEmpty();
        if (made) {
            int m = 0;
            const std::string call = formula_poly_random(
                lo, hi, low, high, QRandomGenerator::global()->generate64(),
                &m, across);
            formula->setText(QString::fromStdString(call));
            formula->setCursorPosition(0);
            status->setText(
                QObject::tr("Length %1, degree %2.").arg(m).arg(m - 1));
        } else {
            formula->clear();
            status->setText(why);
        }
        again->setEnabled(made);
        copy->setEnabled(made);
    };
    for (QLineEdit *field : {shortest, longest, smallest, largest, delta})
        QObject::connect(field, &QLineEdit::textChanged, box, generate);
    QObject::connect(again, &QPushButton::clicked, box, generate);
    QObject::connect(copy, &QPushButton::clicked, box, [formula, status]() {
        QGuiApplication::clipboard()->setText(formula->text());
        status->setText(QObject::tr("Copied to the clipboard: paste it into "
                                    "Fractal -> User formula."));
    });
    generate();

    QWidget *page = new QWidget;
    QVBoxLayout *layout = new QVBoxLayout(page);
    layout->addWidget(box);
    layout->addStretch();
    return page;
}

void ui_formulahelp(struct uih_context * /*uih*/)
{
    ui_formulahelp_open(-1);
}

void ui_formulahelp_open(int tab)
{
    static QWidget *window = NULL;
    static QTabWidget *tabs = NULL;

    if (window == NULL) {
        window = new QWidget(NULL, Qt::Window);
        window->setWindowTitle(
            QObject::tr("%1 - user formula reference").arg(XaoS_NAME));
        window->setAttribute(Qt::WA_QuitOnClose, false);
        window->resize(720, 560);

        tabs = new QTabWidget(window);
        tabs->addTab(buildTable(formula_help_functions, QObject::tr("Function"),
                                true),
                     QObject::tr("Functions"));
        tabs->addTab(buildTable(formula_help_variables, QObject::tr("Name")),
                     QObject::tr("Variables"));
        tabs->addTab(buildTable(formula_help_notation, QObject::tr("Written as")),
                     QObject::tr("Notation"));
        tabs->addTab(buildTable(formula_help_values, QObject::tr("Value")),
                     QObject::tr("Values"));
        tabs->addTab(buildTilings(), QObject::tr("Tilings"));
        tabs->addTab(buildUtilities(), QObject::tr("Utilities"));

        QVBoxLayout *layout = new QVBoxLayout(window);
        layout->addWidget(tabs);
    }

    if (tab >= 0 && tab < tabs->count())
        tabs->setCurrentIndex(tab);
    window->show();
    window->raise();
    window->activateWindow();
}
