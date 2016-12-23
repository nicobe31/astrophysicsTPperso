number = 10;
name = num2str(number,'im/IMG_%.4d.CR2');

[subImage] = imBox(name);

A = angle(subImage);