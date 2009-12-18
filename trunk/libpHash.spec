Name:           pHash
Version:        0.7.0
Release:        1%{?dist}
Summary:        pHash, the open source perceptual hashing library

Group:          System Environment/Libraries
License:        GPLv3
URL:            http://www.phash.org/
Source0:        http://www.phash.org/releases/%{name}-%{version}.tar.gz
BuildRoot:      %{_tmppath}/%{name}-%{version}-%{release}-root-%(%{__id_u} -n)

BuildRequires:  fftw-devel >= 3
BuildRequires:  ffmpeg-devel >= 0.5
BuildRequires:  libjpeg-devel
BuildRequires:  libpng-devel
Requires:	fftw >= 3
Requires:	ffmpeg >= 0.5
Requires:	libjpeg
Requires:	libpng

%description
pHash is a perceptual hashing library that allows you to find similar 
media files without them having to be bit-for-bit identical.

%package        devel
Summary:        Development files for %{name}
Group:          Development/Libraries
Requires:       %{name} = %{version}-%{release}

%description    devel
The %{name}-devel package contains libraries and header files for
developing applications that use %{name}.


%prep
%setup -q


%build
%configure --disable-static
make %{?_smp_mflags}


%install
rm -rf $RPM_BUILD_ROOT
make install DESTDIR=$RPM_BUILD_ROOT
find $RPM_BUILD_ROOT -name '*.la' -exec rm -f {} ';'


%clean
rm -rf $RPM_BUILD_ROOT


%post -p /sbin/ldconfig

%postun -p /sbin/ldconfig


%files
%defattr(-,root,root,-)
%doc
%{_libdir}/*.so.*
/usr/bin/add_mvptree
/usr/bin/add_mvptree_audio
/usr/bin/build_mvptree
/usr/bin/build_mvptree_audio
/usr/bin/query_mvptree
/usr/bin/query_mvptree_audio
/usr/bin/test_audio
/usr/bin/test_image
/usr/bin/test_video
/usr/lib/pkgconfig/pHash.pc

%files devel
%defattr(-,root,root,-)
%doc
%{_includedir}/*
%{_libdir}/*.so


%changelog
