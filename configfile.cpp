
#include <iostream>

#include <QString>
#include <QFile>
#include <QTextStream>


#include "configfile.h"

using namespace std;


InputData::InputData(): fname("asxp.conf") {
}

InputData::InputData(QString fnamea): fname(fnamea) {
}

void InputData::set_fname(QString fnamea) {

	fname = fnamea;

}

void InputData::readSettings(int & cnt_blocks) {

	cnt_blocks = 0;

    QFile infile(fname);
    if (!infile.open(QIODevice::ReadOnly | QIODevice::Text)) {
        return;
    }

    QTextStream in(&infile);
    while (!in.atEnd()) {
        QString line = in.readLine();

        if (line.trimmed().length() == 0) {
        	continue;
        }

        if (line.mid(0, 4) == "surf") {
        	QString surftitle = line.mid(5).trimmed();
        	QString surffun = "";

        	QString line1;

        	cout << "surftitle = " << surftitle.toLatin1().data() << endl;

        	// skip empty lines
        	// break on not empty line found
        	while(!in.atEnd()) {
        		line1 = in.readLine();
        		cout << "line1 = " << line1.toLatin1().data() << endl;
        		cout << "line1.trimmed.len" << line1.trimmed().length() << endl;
        		if (line1.trimmed().length() == 0) {
        			continue;
        		} else {
        			break;
        		}
        	}

        	// add first line of poly to surffun
        	surffun += line1.trimmed();

        	// add all following non empty lines to surffun
        	while(!in.atEnd()) {
        		line1 = in.readLine();
        		cout << "line1_1 = " << line1.toLatin1().data() << endl;
        		if (line1.trimmed().length() == 0) {
        			break;
        		}
        		surffun += line1.trimmed();
        	}
        	SettingsBlock surfsettings;

        	cout << "surffun = " << surffun.toLatin1().data() << endl;
        	surfsettings.block_name = surftitle;
        	surfsettings.smap["poly"] = surffun;

        	while(!in.atEnd()) {
        		line1 = in.readLine();
        		if (line1.trimmed().length() > 0) {
        			if (line1.mid(0,8) == "manifold") {
        				QString manifold_sel = line1.mid(9).trimmed();
        				surfsettings.smap["manifold"] = manifold_sel;
        				cout << "manifold = " << manifold_sel.toLatin1().data() << endl;
        			}
        		} else {
        			break;
        		}
        	}

        	settings_list.push_back(surfsettings);

        	++cnt_blocks;

        }
    }
}


void InputData::getSurfBlocks(QList<QString>& surf_block_list) {

	for(SettingsList::iterator it = settings_list.begin(); it != settings_list.end(); ++it) {
		surf_block_list.push_back(it->block_name);
	}

}

void InputData::getBlock(QString blockname, QString& poly) {

	poly = "";

	for(SettingsList::iterator it = settings_list.begin(); it != settings_list.end(); ++it) {
		if (it->block_name == blockname) {
			poly = it->smap["poly"];
		}
	}
}


void InputData::reset()
{
	settings_list.clear();
}
