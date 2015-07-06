<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Extract a isocontour</description>
    <executable>FieldConvert</executable>
    <parameters>-m isocontour:fieldstr="u+v":fieldvalue=0.5:fieldname="UplusV":smooth chan3D.xml chan3D.fld isocontour.dat</parameters>
    <files>
        <file description="Session File">chan3D.xml</file>
        <file description="Session File">chan3D.fld</file>
    </files>
    <metrics>
        <metric type="file" id="1">
            <file filename="isocontour.dat">
                <sha1>90dcaca2feb7d4af87d715b90cf4097ecd5e82eb</sha1>
             </file>
         </metric>
    </metrics>
</test>

