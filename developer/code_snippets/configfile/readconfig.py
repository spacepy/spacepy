## this is an example configuration reading tool
## MOdify and add to whatever as needed
## Brian Larsen 24-May-2010

import ConfigParser
def readconfig(my_cfg, config_filepath):
    # Create a ConfigParser object, to read the config file
    cfg=ConfigParser.ConfigParser()
    fp=open(config_filepath, "rb")
    cfg.readfp(fp)
    fp.close()
    # Read each parameter in turn
    # Web section values
    if cfg.has_option("web", "port"):
        my_cfg["PORT"]=cfg.get("web", "port")
    if cfg.has_option("web", "browser_refresh"):
        my_cfg["BROWSER_REFRESH"]=cfg.get("web", "browser_refresh")
    # service section values
    if cfg.has_option("service", "service_name"):
        my_cfg["SERVICE_NAME"]=cfg.get("service", "service_name")
    if cfg.has_option("service", "poll_interval"):
        my_cfg["POLL_INTERVAL"]=cfg.get("service", "poll_interval")
    if cfg.has_option("service", "log_records"):
        my_cfg["LOG_RECORDS"]=cfg.get("service", "log_records")
    return

# usage
# a = {}
# readconfig(a, 'test.cfg')


