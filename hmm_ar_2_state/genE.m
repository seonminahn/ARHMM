function a = genE(nvars)

a = rand(nvars)+4*eye(nvars);
for i=1:nvars
    a(i,:)=a(i,:)/sum(a(i,:));
end