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

/* One tab: a two-column table of names and descriptions, with the section
 * rows spanning both columns.
 *
 * QTableWidget rather than anything lighter because it scrolls by itself and
 * keeps the header in place while it does, which is what makes a list of
 * ninety-six entries usable. */
static QWidget *buildTable(const struct formula_help_row *rows,
                           const QString &nameHeading)
{
    int count = 0;
    for (const struct formula_help_row *r = rows; r->name || r->section; r++)
        count++;

    QTableWidget *table = new QTableWidget(count, 2);
    table->setHorizontalHeaderLabels(QStringList()
                                     << nameHeading
                                     << QObject::tr("What it does"));
    table->verticalHeader()->setVisible(false);
    table->setEditTriggers(QAbstractItemView::NoEditTriggers);
    table->setSelectionBehavior(QAbstractItemView::SelectRows);
    table->setAlternatingRowColors(true);
    table->setWordWrap(true);
    table->horizontalHeader()->setSectionResizeMode(0,
                                                    QHeaderView::ResizeToContents);
    table->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Stretch);

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
            /* The heading is one cell across both columns. */
            table->setSpan(row, 0, 1, 2);
            continue;
        }
        QTableWidgetItem *name = new QTableWidgetItem(QString(r->name));
        QFont fixed = QFontDatabase::systemFont(QFontDatabase::FixedFont);
        name->setFont(fixed);
        table->setItem(row, 0, name);
        table->setItem(row, 1, new QTableWidgetItem(QObject::tr(r->summary)));
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
        tabs->addTab(buildTable(formula_help_functions, QObject::tr("Function")),
                     QObject::tr("Functions"));
        tabs->addTab(buildTable(formula_help_variables, QObject::tr("Name")),
                     QObject::tr("Variables"));
        tabs->addTab(buildTable(formula_help_notation, QObject::tr("Written as")),
                     QObject::tr("Notation"));

        QVBoxLayout *layout = new QVBoxLayout(window);
        layout->addWidget(tabs);
    }

    window->show();
    window->raise();
    window->activateWindow();
}
