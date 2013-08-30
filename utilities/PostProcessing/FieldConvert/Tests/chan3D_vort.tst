<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Calculate Vorticity using Fieldconvert </description>
    <executable>Fieldconvert</executable>
    <parameters>-m vorticity chan3D.xml chan3D.fld chan_vort.fld</parameters>
    <files>
        <file description="Session File">chan3D.xml</file>
	<file description="Session File">chan3D.fld</file>
	<file description="Session File">chan3D_vort.fld</file>
    </files>
    <metrics>
      <file filename="chan3D_vort.fld">
        <sha1>sha1 hash here</sha1>
      </file>
    </metrics>
</test>

