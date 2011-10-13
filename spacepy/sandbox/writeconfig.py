## this is an example configuration writing tool
## MOdify and add to whatever as needed
## Brian Larsen 24-May-2010

import ConfigParser
def writeconfig(my_cfg, config_filepath):
    cfg=ConfigParser.ConfigParser()
    cfg.add_section("web")
    cfg.add_section("service")
    cfg.set("web", "port", my_cfg["PORT"])
    cfg.set("web", "browser_refresh", my_cfg["BROWSER_REFRESH"])
    cfg.set("service", "service_name", my_cfg["SERVICE_NAME"])
    cfg.set("service", "poll_interval", my_cfg["POLL_INTERVAL"])
    cfg.set("service", "log_records", my_cfg["LOG_RECORDS"])
    fp=open(config_filepath, "wb")
    cfg.write(fp)
    fp.close()
    return

# usage
# writeconfig({'PORT':12, 'BROWSER_REFRESH':False, 'SERVICE_NAME':'daman', 'POLL_INTERVAL':1200, 'LOG_RECORDS':True}, 'test.cfg')