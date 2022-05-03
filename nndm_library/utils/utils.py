def get_public_attr(obj):
    # obtain the attr to be explored
    attr = dir(obj)

    # filtrate the ones that are public
    attr = [atr for atr in attr if atr[:2] != '__' and atr[-2:] != '__']

    return attr