<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.1//EN" "http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd"><html xmlns="http://www.w3.org/1999/xhtml"><!-- InstanceBegin template="/Templates/Main.dwt" codeOutsideHTMLIsLocked="false" -->
  <head>
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
    <link rel="shortcut icon" href="images/favicon.ico" >
    <!-- InstanceBeginEditable name="doctitle" -->
    <title>Nektar++ - Downloads</title>
    <!-- InstanceEndEditable -->
 <link rel="stylesheet" href="style.css">
    <!-- InstanceBeginEditable name="CSS_Style" -->
    <style type="text/css">
    <!--
    #menu_Downloads {
      font-weight: bold;
    }
    -->
    </style>
    <!-- InstanceEndEditable -->
    <!--[if IE]>
    <style type="text/css"> 
      .twoColHybLt #sidebar1 { padding-top: 30px; margin-top: 20px; }
      .twoColHybLt #mainContent { zoom: 1; padding-top: 15px; }
    </style>
    <![endif]-->
    <!-- InstanceParam name="LibraryMenu" type="boolean" value="false" -->
    <!-- InstanceParam name="CompileMenu" type="boolean" value="false" -->
  </head>
  <body class="twoColHybLt">
    <div id="container">
      <div id="sidebar1">
        <h3 style="text-align: center;"><a href="index.html" id="menu_Nektar"><img src="images/nektar.png" alt="Nektar++" width="155" height="38" /></a></h3>

<p><a href="downloads.html" id="menu_Downloads">Downloads</a></p>
<p><a href="compile.html" id="menu_Compile">Compilation Instructions</a></p>
<p><a href="usage.html" id="menu_Usage">Example Usage</a></p>
<p><a href="educational_material.html" id="menu_EducationalMaterial">Educational
Material</a></p>
<p><a href="code/index.html" id="menu_Library">Documentation</a></p>
<p><a href="team.html" id="menu_Team">Team Members</a></p>
<p><a href="publications.html" id="menu_Publications">Publications</a></p>
<p><a href="license.html" id="menu_License">License</a></p>
<p><a href="acknowledgments.html" id="menu_Acknowledgments">Acknowledgments</   a></p>
<p><a href="mailto:nektar-inquiry@sci.utah.edu" id="menu_Contact">Contact</a></ p>
      </div>
      <div id="mainContent">
        <!-- InstanceBeginEditable name="Title" -->
        <h1>Downloads</h1>
    	<!-- InstanceEndEditable -->
        <hr>
	<ul>
      <li>
         <h3 align="left">release 3.1.0 (January 2012) </h3>
         <p><a href="downloads/nektar++-3.1.0.tar.gz">Nektar++ Linux/Windows Code</a></p>
         <p><a href="downloads/ThirdParty-3.1.0.tar.gz">Required Third Party Libraries</a></p>
      </li>
      <li>
         <h3 align="left">release 3.0.1 (October 2011) </h3>
         <p><a href="downloads/nektar++-3.0.1.tar.gz">Nektar++ Linux/Windows Code</a></p>
         <p><a href="downloads/ThirdParty-3.0.1.tar.gz">Required Third Party Libraries</a></p>
      </li>
      <li>
         <h3 align="left">release 3.0.0 (August 2011) </h3>
         <p><a href="downloads/nektar++-3.0.0.tar.gz">Nektar++ Linux/Windows Code</a></p>
         <p><a href="downloads/ThirdParty-3.0.0.tar.gz">Required Third Party Libraries</a></p>
      </li>
      <br/>
      <li>
         <h3 align="left">release 2.0 (October 2010) </h3>
         <p><a href="downloads/Nektar++-2.0.tar.gz">Nektar++ Linux/Windows Code</a></p>
         <p><a href="downloads/ThirdParty-2.0.tar.gz">Required Third Party Libraries</a></p>
      </li>
      <br/>
	  <li>
	    <h3 align="left">release 1.1 (May 2009) </h3>
	    <p><a href="downloads/Nektar++1.1.tar.gz">Nektar++ Linux/Windows Code</a></p>
    	    <p><a href="downloads/ThirdParty1.1.tar.gz">Required Third Party Libraries</a></p>
	  </li>
	  <BR>
	  <li>
	    <h3 align="left">release 1.0 (May 2008) </h3>
	    <p><a href="downloads/Nektar++1.0.tar.gz">Nektar++ Linux/Windows Code</a></p>
    	    <p><a href="downloads/ThirdParty1.0.tar.gz">Required Third Party Libraries</a></p>
	  </li>
	</ul>

      </div>
    </div>
  </body>
<!-- InstanceEnd --></html>
<!--
<?php
// Load connection constants.
require('common.php');

// Connect to MySQL.
$link = mysql_connect($hostname, $username, $password);
if (!$link) exit(1);

// Select the correct database.
$db_selected = mysql_select_db($database, $link);
if (!$db_selected) exit(1);

// Add the person to the database.
$query = sprintf("INSERT INTO `$table` (`Name`, `Email`, `Affiliation`) VALUES ('%s', '%s', '%s')",
                 mysql_real_escape_string($_GET["Name"], $link),
                 mysql_real_escape_string($_GET["Email"], $link),
                 mysql_real_escape_string($_GET["Affiliation"], $link));
mysql_query($query, $link);

// Close the connection to MySQL.
mysql_close($link);
?>
-->
