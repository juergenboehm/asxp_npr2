#ifndef __configfile_h
#define __configfile_h

#include <QMap>
#include <QList>
#include <QString>


class InputData {

public:
	InputData();
	InputData(QString fnamea);

	void set_fname(QString fnamea);

	void readSettings(int & cnt_blocks);

	void getSurfBlocks(QList<QString> & surf_block_list);

	void getBlock(QString blockname, QString & poly);

	void reset();

	QString fname;

	typedef QMap<QString, QString> SettingsMap;

	typedef struct {
		QString block_name;
		SettingsMap smap;
	} SettingsBlock;

	typedef QList<SettingsBlock> SettingsList;

	SettingsList settings_list;

};


#endif
