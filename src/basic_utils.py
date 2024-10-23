import os

# make dir
def make_dir_check(dirPath):
    """
    Ensures that a directory exists at the specified path. If the directory does not exist, it is created.

    Args:
        dirPath (str): The path of the directory to check or create.

    Returns:
        None

    Side Effects:
        - Creates the directory at the specified path if it does not exist.
        - Prints a message indicating whether the directory was created or already existed.
    """
    dirPath = dirPath.strip()
    dirPath = dirPath.rstrip("\\")
    isExist = os.path.exists(dirPath)
    if not isExist:
        os.makedirs(dirPath, exist_ok=True)
    else:
        print(f"{dirPath} already exists")
    print(f"Directory created: {dirPath}")