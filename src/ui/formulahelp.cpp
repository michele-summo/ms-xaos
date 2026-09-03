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
        table->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Fixed);
        table->setColumnWidth(1,
                              taken + 3 * metrics.horizontalAdvance(QChar('m')));
    }
    table->horizontalHeader()->setSectionResizeMode(columns - 1,
                                                    QHeaderView::Stretch);

    int row = 0;
    for (const struct formula_help_row *r = rows; r->name || r->section;
         r++, row++) {
        if (r->section != NULL) {
            QTableWidgetItem *item =
                new QTableWidgetItem(QObject::tr(r->section));
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
                       new QTableWidgetItem(QObject::tr(r->summary)));
    }
    table->resizeRowsToContents();
    return table;
}

void ui_formulahelp(struct uih_context * /*uih*/)
{
    static QWidget *window = NULL;

    if (window == NULL) {
        window = new QWidget(NULL, Qt::Window);
        window->setWindowTitle(
            QObject::tr("%1 - user formula reference").arg(XaoS_NAME));
        window->setAttribute(Qt::WA_QuitOnClose, false);
        window->resize(720, 560);

        QTabWidget *tabs = new QTabWidget(window);
        tabs->addTab(buildTable(formula_help_functions, QObject::tr("Function"),
                                true),
                     QObject::tr("Functions"));
        tabs->addTab(buildTable(formula_help_variables, QObject::tr("Name")),
                     QObject::tr("Variables"));
        tabs->addTab(buildTable(formula_help_notation, QObject::tr("Written as")),
                     QObject::tr("Notation"));
        tabs->addTab(buildTable(formula_help_values, QObject::tr("Value")),
                     QObject::tr("Values"));

        QVBoxLayout *layout = new QVBoxLayout(window);
        layout->addWidget(tabs);
    }

    window->show();
    window->raise();
    window->activateWindow();
}
