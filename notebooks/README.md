如果你使用新的从full image crop的图作为not digit的train data, 请使用commit：second version.
否则使用其他commit。

使用full image crop做训练的好处是，在全图detection的时候，在与数字没有交集的地方出现无检测的可能性下降很多（相比于cifar）。

