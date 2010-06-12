Name:           pHash
Version:        0.9.1
Release:        1%{?dist}
Summary:        pHash, the open source perceptual hashing library

Group:          System Environment/Libraries
License:        GPLv3
URL:            http://www.phash.org/
Source0:        http://www.phash.org/releases/%{name}-%{version}.tar.gz
BuildRoot:      %{_tmppath}/%{name}-%{version}-%{release}-root-%(%{__id_u} -n)

BuildRequires:  ffmpeg-devel >= 0.5
BuildRequires:  libjpeg-devel
BuildRequires:  libpng-devel
BuildRequires:  libsndfile-devel
BuildRequires:  libsamplerate-devel
Requires:	ffmpeg >= 0.5
Requires:	libjpeg
Requires:	libpng
Requires:  	libsndfile
Requires:  	libsamplerate

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
/usr/lib/pkgconfig/pHash.pc

%files devel
%defattr(-,root,root,-)
%doc
%{_includedir}/*
%{_libdir}/*.so


%changelog
