def flatten_dict(dictionary):
    """
    Recursively flattens a nested dictionary. It only retains the keys at their
    deepest point in the dictionary.
    """
    nested_check = is_dict_nested(dictionary)

    if nested_check:
        result = {}
        for key, value in dictionary.items():
            if isinstance(value, dict):
                # If value is a dictionary, recursively flatten it
                # and update our result dictionary
                nested_flat = flatten_dict(value)
                result.update(nested_flat)
            else:
                # Add the key-value pair to the result
                result[key] = value
        return result
    else:
        return dictionary


def is_dict_nested(dictionary):
    """
    Check if a dictionary is nested (contains at least one dictionary as a
    value).
    """
    if not isinstance(dictionary, dict):
        return False
    
    for value in dictionary.values():
        if isinstance(value, dict):
            return True
    return False