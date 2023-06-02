from .mgw import mgw
from .prot import prot
from .helt import helt
from .roll import roll
import numpy as np
models = {"mgw": mgw.MGW, "prot": prot.PROT, "helt": helt.HELT, "roll": roll.ROLL}
for feat in models.keys():
    for length in models[feat].keys():
        models[feat][length]["frequencies"] = np.array(models[feat][length]["frequencies"], dtype=np.dtype('d'))
        models[feat][length]["bins"] = np.array(models[feat][length]["bins"], dtype=np.dtype('d'))
    