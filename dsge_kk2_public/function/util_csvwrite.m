function [] = util_csvwrite(filename, x, str)

fid = fopen(filename,'wt'); % �������ݗp�Ƀt�@�C���I�[�v��
[rows,cols] = size(x);
if cols > 1
  fprintf(fid, '%s,', str{1,1:end-1}); % ������̏����o��
  fprintf(fid, '%s\n', str{1,end}); % �s���̕�����́A���s���܂߂ďo��
  for i = 1:rows
      fprintf(fid, '%f,', x(i,1:end-1));
      fprintf(fid, '%f\n', x(i, end));
  end
else
  fprintf(fid, '%s\n', str{1});
  for i = 1:rows
      fprintf(fid, '%f\n', x(i));
  end
end
fclose('all');
 % �t�@�C���N���[�Y
